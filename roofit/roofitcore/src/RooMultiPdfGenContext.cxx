/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * Authors:                                                                  *
 *   Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooMultiPdfGenContext.cxx
\class RooMultiPdfGenContext
\ingroup Roofitcore

Generator context for a RooMultiPdf that is used when the index category of
the RooMultiPdf is among the observables to be generated. In that case, the
RooMultiPdf is treated as a mixture: for each event the index category is
sampled according to the relative expected yields of the component pdfs, and
the observables are generated from the selected component. This mirrors the
behaviour of RooSimGenContext for RooSimultaneous.

If no prototype data is given, the relative fractions are computed from the
extended term (expected events) of the component pdfs, which therefore need to
be extended.
**/

#include "RooMultiPdfGenContext.h"

#include "RooMultiPdf.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooGlobalFunc.h"

#include <string>

////////////////////////////////////////////////////////////////////////////////
/// Constructor of a generator context for a RooMultiPdf. A dedicated generator
/// context is created for each component pdf and generation of events is
/// delegated to the appropriate component context.

RooMultiPdfGenContext::RooMultiPdfGenContext(const RooMultiPdf &model, const RooArgSet &vars,
                                             const RooDataSet *prototype, const RooArgSet *auxProto, bool verbose)
   : RooAbsGenContext(model, vars, prototype, auxProto, verbose), _pdf(&model)
{
   auto const &idxCat = static_cast<RooAbsCategoryLValue const &>(model.indexCategory().arg());

   RooArgSet pdfVars(vars);
   RooArgSet allPdfVars(pdfVars);
   if (prototype)
      allPdfVars.add(*prototype->get(), true);

   // The index category must be among the observables to generate
   if (!allPdfVars.find(idxCat.GetName())) {
      oocoutE(_pdf, Generation) << "RooMultiPdfGenContext::ctor(" << GetName()
                                << ") ERROR: the index category must be among the observables to generate"
                                << std::endl;
      _isValid = false;
      _numPdf = 0;
      _haveIdxProto = false;
      return;
   }

   // We must either have the prototype or an extended likelihood to determine
   // the relative fractions of the components
   _haveIdxProto = prototype != nullptr;
   _idxCatName = idxCat.GetName();

   _numPdf = model.getNumPdfs();

   if (!_haveIdxProto) {
      for (int i = 0; i < _numPdf; ++i) {
         if (!model.getPdf(i)->canBeExtended()) {
            oocoutE(_pdf, Generation)
               << "RooMultiPdfGenContext::ctor(" << GetName() << ") ERROR: Need either extended component pdfs"
               << " or prototype data to calculate number of events per component" << std::endl;
            _isValid = false;
            _numPdf = 0;
            return;
         }
      }
   }

   // Initialize fraction threshold array (used only in extended mode)
   _fracThresh.resize(_numPdf + 1);
   _fracThresh[0] = 0;

   // Create a generator context for each component pdf
   _allVarsPdf.add(allPdfVars);
   for (int i = 0; i < _numPdf; ++i) {
      RooAbsPdf *pdf = model.getPdf(i);

      _gcList.emplace_back(pdf->genContext(pdfVars, prototype, auxProto, verbose));
      _gcList.back()->SetName(pdf->GetName());

      // In a RooMultiPdf the component at list position i is selected by the
      // index category value i (see RooMultiPdf::getPdf).
      _gcIndex.push_back(i);

      _fracThresh[i + 1] = _fracThresh[i] + (_haveIdxProto ? 0 : pdf->expectedEvents(&allPdfVars));
   }

   // Normalize fraction threshold array
   if (!_haveIdxProto) {
      for (int i = 0; i < _numPdf; ++i)
         _fracThresh[i] /= _fracThresh[_numPdf];
   }

   // Clone the index category
   _idxCatSet = std::make_unique<RooArgSet>();
   RooArgSet(idxCat).snapshot(*_idxCatSet, true);
   if (_idxCatSet->empty()) {
      oocoutE(_pdf, Generation) << "RooMultiPdfGenContext::RooMultiPdfGenContext(" << GetName()
                                << ") Couldn't deep-clone index category, abort," << std::endl;
      throw std::string("RooMultiPdfGenContext::RooMultiPdfGenContext() Couldn't deep-clone index category, abort");
   }

   _idxCat = static_cast<RooAbsCategoryLValue *>(_idxCatSet->find(idxCat.GetName()));
}

////////////////////////////////////////////////////////////////////////////////
/// Attach the index category clone to the given event buffer

void RooMultiPdfGenContext::attach(const RooArgSet &args)
{
   if (_idxCat->isDerived()) {
      _idxCat->recursiveRedirectServers(args);
   }

   // Forward attach call to all components
   for (auto &item : _gcList) {
      item->attach(args);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform one-time initialization of the generator context

void RooMultiPdfGenContext::initGenerator(const RooArgSet &theEvent)
{
   // Attach the index category clone to the event
   if (_idxCat->isDerived()) {
      _idxCat->recursiveRedirectServers(theEvent);
   } else {
      _idxCat = static_cast<RooAbsCategoryLValue *>(theEvent.find(_idxCat->GetName()));
   }

   // Update fractions reflecting possible new parameter values
   updateFractions();

   // Forward initGenerator call to all components
   for (auto &item : _gcList) {
      item->initGenerator(theEvent);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create an empty dataset to hold the events that will be generated

RooDataSet *RooMultiPdfGenContext::createDataSet(const char *name, const char *title, const RooArgSet &obs)
{
   // If the observables do not contain the index, make a plain dataset
   if (!obs.contains(*_idxCat)) {
      return new RooDataSet(name, title, obs);
   }

   if (!_protoData) {
      std::map<std::string, RooAbsData *> dmap;
      for (int i = 0; i < _numPdf; ++i) {
         std::string label = _idxCat->lookupName(_gcIndex[i]);
         RooAbsPdf *slicePdf = _pdf->getPdf(i);
         std::unique_ptr<RooArgSet> sliceObs{slicePdf->getObservables(obs)};
         std::string sliceName = name + ("_slice_" + label);
         std::string sliceTitle = title + (" (index slice " + label + ")");
         dmap[label] = new RooDataSet(sliceName, sliceTitle, *sliceObs);
      }
      using namespace RooFit;
      _protoData = std::make_unique<RooDataSet>(name, title, obs, Index(static_cast<RooCategory &>(*_idxCat)),
                                                Link(dmap), OwnLinked());
   }

   return new RooDataSet(*_protoData, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an event. The index state is either taken from the prototype data
/// or sampled from the fraction threshold table, and the observables are
/// generated by the corresponding component generator context.

void RooMultiPdfGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
   if (_haveIdxProto) {

      // Lookup pdf from selected prototype index state
      Int_t gidx(0);
      Int_t cidx = _idxCat->getCurrentIndex();
      for (Int_t i = 0; i < (Int_t)_gcIndex.size(); i++) {
         if (_gcIndex[i] == cidx) {
            gidx = i;
            break;
         }
      }
      RooAbsGenContext *cx = _gcList[gidx].get();
      if (cx) {
         cx->generateEvent(theEvent, remaining);
      } else {
         oocoutW(_pdf, Generation) << "RooMultiPdfGenContext::generateEvent: WARNING, no PDF to generate event of type "
                                   << cidx << std::endl;
      }

   } else {

      // Throw a random number and select the component from the fraction
      // threshold table, i.e. pick the first component whose upper threshold
      // exceeds the random number (inversion of the cumulative distribution).
      double rand = RooRandom::uniform();
      for (Int_t i = 0; i < _numPdf; i++) {
         if (rand < _fracThresh[i + 1]) {
            _gcList[i]->generateEvent(theEvent, remaining);
            // Write through to the index category because it is written to the dataset
            _idxCat->setIndex(_gcIndex[i]);
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute the fraction threshold array to reflect the current parameter values

void RooMultiPdfGenContext::updateFractions()
{
   if (_haveIdxProto)
      return;

   for (int i = 0; i < _numPdf; ++i) {
      _fracThresh[i + 1] = _fracThresh[i] + _pdf->getPdf(i)->expectedEvents(&_allVarsPdf);
   }

   // Normalize fraction threshold array
   for (int i = 0; i < _numPdf; ++i)
      _fracThresh[i] /= _fracThresh[_numPdf];
}

////////////////////////////////////////////////////////////////////////////////
/// Set the traversal order of the prototype data. This information is passed to
/// all component generator contexts.

void RooMultiPdfGenContext::setProtoDataOrder(Int_t *lut)
{
   RooAbsGenContext::setProtoDataOrder(lut);

   for (auto &item : _gcList) {
      item->setProtoDataOrder(lut);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooMultiPdfGenContext::printMultiline(std::ostream &os, Int_t content, bool verbose, TString indent) const
{
   RooAbsGenContext::printMultiline(os, content, verbose, indent);
   os << indent << "--- RooMultiPdfGenContext ---" << std::endl;
   os << indent << "Using PDF ";
   _pdf->printStream(os, kName | kArgs | kClassName, kSingleLine, indent);
   os << indent << "List of component generators" << std::endl;

   TString indent2(indent);
   indent2.Append("    ");

   for (auto &item : _gcList) {
      item->printMultiline(os, content, verbose, indent2);
   }
}
