/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooSimultaneous.cxx
\class RooSimultaneous
\ingroup Roofitcore

Facilitates simultaneous fitting of multiple PDFs to subsets of a given dataset.
The class takes an index category, which is used as a selector
for PDFs, and a list of PDFs, each associated
with a state of the index category. RooSimultaneous always returns
the value of the PDF that is associated with the current value
of the index category.

Extended likelihood fitting is supported if all components support
extended likelihood mode. The expected number of events by a RooSimultaneous
is that of the component p.d.f. selected by the index category.

The index category can be accessed using indexCategory().

###Generating events
When generating events from a RooSimultaneous, the index category has to be added to
the dataset. Further, the PDF needs to know the relative probabilities of each category, i.e.,
how many events are in which category. This can be achieved in two ways:
- Generating with proto data that have category entries: An event from the same category as
in the proto data is created for each event in the proto data.
See RooAbsPdf::generate(const RooArgSet&,const RooDataSet&,Int_t,bool,bool,bool) const.
- No proto data: A category is chosen randomly.
\note This requires that the PDFs building the simultaneous are extended. In this way,
the relative probability of each category can be calculated from the number of events
in each category.
**/

#include "RooSimultaneous.h"

#include "Roo1DTable.h"
#include "RooAbsCategoryLValue.h"
#include "RooAbsData.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinSamplingPdf.h"
#include "RooBinWidthFunction.h"
#include "RooCategory.h"
#include "RooCmdConfig.h"
#include "RooCompositeDataStore.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooNameReg.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimGenContext.h"
#include "RooSimSplitGenContext.h"
#include "RooSuperCategory.h"

#include "RooFitImplHelpers.h"

#include <RooFit/Detail/RooChannelIndicatorPdf.h>

#include <ROOT/StringUtils.hxx>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string_view>

namespace {

std::map<std::string, RooAbsPdf *> createPdfMap(const RooArgList &inPdfList, RooAbsCategoryLValue &inIndexCat)
{
   std::map<std::string, RooAbsPdf *> pdfMap;
   auto indexCatIt = inIndexCat.begin();
   for (unsigned int i = 0; i < inPdfList.size(); ++i) {
      auto pdf = static_cast<RooAbsPdf *>(&inPdfList[i]);
      const auto &nameIdx = (*indexCatIt++);
      pdfMap[nameIdx.first] = pdf;
   }
   return pdfMap;
}

void replaceOrAdd(RooLinkedList &lst, TObject &obj)
{
   TObject *old = lst.FindObject(obj.GetName());
   if (old)
      lst.Replace(old, &obj);
   else
      lst.Add(&obj);
}

} // namespace

RooSimultaneous::InitializationOutput::~InitializationOutput() = default;

void RooSimultaneous::InitializationOutput::addPdf(const RooAbsPdf &pdf, std::string const &catLabel)
{
   finalPdfs.push_back(&pdf);
   finalCatLabels.emplace_back(catLabel);
}

using std::string;




////////////////////////////////////////////////////////////////////////////////
/// Constructor with index category. PDFs associated with indexCat
/// states can be added after construction with the addPdf() function.
///
/// RooSimultaneous can function without having a PDF associated
/// with every single state. The normalization in such cases is taken
/// from the number of registered PDFs, but getVal() will assert if
/// when called for an unregistered index state.

RooSimultaneous::RooSimultaneous(const char *name, const char *title,
             RooAbsCategoryLValue& inIndexCat) :
  RooSimultaneous{name, title, std::map<std::string, RooAbsPdf*>{}, inIndexCat}
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor from index category and full list of PDFs.
/// In this constructor form, a PDF must be supplied for each indexCat state
/// to avoid ambiguities. The PDFs are associated with the states of the
/// index category as they appear when iterating through the category states
/// with RooAbsCategory::begin() and RooAbsCategory::end(). This usually means
/// they are associated by ascending index numbers.
///
/// PDFs may not overlap (i.e. share any variables) with the index category (function)

RooSimultaneous::RooSimultaneous(const char *name, const char *title,
             const RooArgList& inPdfList, RooAbsCategoryLValue& inIndexCat) :
  RooSimultaneous{name, title, createPdfMap(inPdfList, inIndexCat), inIndexCat}
{
  if (inPdfList.size() != inIndexCat.size()) {
    std::stringstream errMsg;
    errMsg << "RooSimultaneous::ctor(" << GetName()
           << " ERROR: Number PDF list entries must match number of index category states, no PDFs added";
    coutE(InputArguments) << errMsg.str() << std::endl;
    throw std::invalid_argument(errMsg.str());
  }
}


////////////////////////////////////////////////////////////////////////////////

RooSimultaneous::RooSimultaneous(const char *name, const char *title, std::map<string, RooAbsPdf *> pdfMap,
                                 RooAbsCategoryLValue &inIndexCat)
   : RooSimultaneous(name, title, std::move(*initialize(name ? name : "", inIndexCat, pdfMap)))
{
}

/// For internal use in RooFit.
RooSimultaneous::RooSimultaneous(const char *name, const char *title,
                                 RooFit::Detail::FlatMap<std::string, RooAbsPdf *> const &pdfMap,
                                 RooAbsCategoryLValue &inIndexCat)
   : RooSimultaneous(name, title, RooFit::Detail::flatMapToStdMap(pdfMap), inIndexCat)
{
}

RooSimultaneous::RooSimultaneous(const char *name, const char *title, RooSimultaneous::InitializationOutput &&initInfo)
   : RooAbsPdf(name, title),
     _plotCoefNormSet("!plotCoefNormSet", "plotCoefNormSet", this, false, false),
     _partIntMgr(this, 10),
     _indexCat("indexCat", "Index category", this, *initInfo.indexCat)
{
   for (std::size_t i = 0; i < initInfo.finalPdfs.size(); ++i) {
      addPdf(*initInfo.finalPdfs[i], initInfo.finalCatLabels[i].c_str());
   }

   // Take ownership of eventual super category
   if (initInfo.superIndex) {
      addOwnedComponents(std::move(initInfo.superIndex));
   }
}

/// \cond ROOFIT_INTERNAL

// This class cannot be locally defined in initialize as it cannot be
// used as a template argument in that case
namespace RooSimultaneousAux {
  struct CompInfo {
    RooAbsPdf* pdf ;
    RooSimultaneous* simPdf ;
    const RooAbsCategoryLValue* subIndex ;
    std::unique_ptr<RooArgSet> subIndexComps;
  } ;
}

/// \endcond

std::unique_ptr<RooSimultaneous::InitializationOutput>
RooSimultaneous::initialize(std::string const& name, RooAbsCategoryLValue &inIndexCat,
                            std::map<std::string, RooAbsPdf *> const& pdfMap)

{
  auto out = std::make_unique<RooSimultaneous::InitializationOutput>();
  out->indexCat = &inIndexCat;

  // First see if there are any RooSimultaneous input components
  bool simComps(false) ;
  for (auto const& item : pdfMap) {
    if (dynamic_cast<RooSimultaneous*>(item.second)) {
      simComps = true ;
      break ;
    }
  }

  // If there are no simultaneous component p.d.f. do simple processing through addPdf()
  if (!simComps) {
    for (auto const& item : pdfMap) {
      out->addPdf(*item.second,item.first);
    }
    return out;
  }

  std::string msgPrefix = "RooSimultaneous::initialize(" + name + ") ";

  // Issue info message that we are about to do some rearranging
  oocoutI(nullptr, InputArguments) << msgPrefix << "INFO: one or more input component of simultaneous p.d.f.s are"
         << " simultaneous p.d.f.s themselves, rewriting composite expressions as one-level simultaneous p.d.f. in terms of"
         << " final constituents and extended index category" << std::endl;


  RooArgSet allAuxCats ;
  std::map<string,RooSimultaneousAux::CompInfo> compMap ;
  for (auto const& item : pdfMap) {
    RooSimultaneousAux::CompInfo ci ;
    ci.pdf = item.second ;
    RooSimultaneous* simComp = dynamic_cast<RooSimultaneous*>(item.second) ;
    if (simComp) {
      ci.simPdf = simComp ;
      ci.subIndex = &simComp->indexCat() ;
      ci.subIndexComps = simComp->indexCat().isFundamental()
          ? std::make_unique<RooArgSet>(simComp->indexCat())
          : std::unique_ptr<RooArgSet>(simComp->indexCat().getVariables());
      allAuxCats.add(*ci.subIndexComps,true) ;
    } else {
      ci.simPdf = nullptr;
      ci.subIndex = nullptr;
    }
    compMap[item.first] = std::move(ci);
  }

  // Construct the 'superIndex' from the nominal index category and all auxiliary components
  RooArgSet allCats(inIndexCat) ;
  allCats.add(allAuxCats) ;
  std::string siname = name + "_index";
  out->superIndex = std::make_unique<RooSuperCategory>(siname.c_str(),siname.c_str(),allCats) ;
  auto *superIndex = out->superIndex.get();
  out->indexCat = superIndex;

  // Now process each of original pdf/state map entries
  for (auto const& citem : compMap) {

    RooArgSet repliCats(allAuxCats) ;
    if (citem.second.subIndexComps) {
      repliCats.remove(*citem.second.subIndexComps) ;
    }
    inIndexCat.setLabel(citem.first.c_str()) ;

    if (!citem.second.simPdf) {

      // Entry is a plain p.d.f. assign it to every state permutation of the repliCats set
      RooSuperCategory repliSuperCat("tmp","tmp",repliCats) ;

      // Iterator over all states of repliSuperCat
      for (const auto& nameIdx : repliSuperCat) {
        // Set value
        repliSuperCat.setLabel(nameIdx.first) ;
        // Retrieve corresponding label of superIndex
        string superLabel = superIndex->getCurrentLabel() ;
        out->addPdf(*citem.second.pdf,superLabel);
        oocxcoutD(static_cast<RooAbsArg*>(nullptr), InputArguments) << msgPrefix
                << "assigning pdf " << citem.second.pdf->GetName() << " to super label " << superLabel << std::endl ;
      }
    } else {

      // Entry is a simultaneous p.d.f

      if (repliCats.empty()) {

        // Case 1 -- No replication of components of RooSim component are required

        for (const auto& type : *citem.second.subIndex) {
          const_cast<RooAbsCategoryLValue*>(citem.second.subIndex)->setLabel(type.first.c_str());
          string superLabel = superIndex->getCurrentLabel() ;
          RooAbsPdf* compPdf = citem.second.simPdf->getPdf(type.first);
          if (compPdf) {
            out->addPdf(*compPdf,superLabel);
            oocxcoutD(static_cast<RooAbsArg*>(nullptr), InputArguments) << msgPrefix
                    << "assigning pdf " << compPdf->GetName() << "(member of " << citem.second.pdf->GetName()
                    << ") to super label " << superLabel << std::endl ;
          } else {
            oocoutW(nullptr, InputArguments) << msgPrefix << "WARNING: No p.d.f. associated with label "
                << type.second << " for component RooSimultaneous p.d.f " << citem.second.pdf->GetName()
                << "which is associated with master index label " << citem.first << std::endl ;
          }
        }

      } else {

        // Case 2 -- Replication of components of RooSim component are required

        // Make replication supercat
        RooSuperCategory repliSuperCat("tmp","tmp",repliCats) ;

        for (const auto& stype : *citem.second.subIndex) {
          const_cast<RooAbsCategoryLValue*>(citem.second.subIndex)->setLabel(stype.first.c_str());

          for (const auto& nameIdx : repliSuperCat) {
            repliSuperCat.setLabel(nameIdx.first) ;
            const string superLabel = superIndex->getCurrentLabel() ;
            RooAbsPdf* compPdf = citem.second.simPdf->getPdf(stype.first);
            if (compPdf) {
              out->addPdf(*compPdf,superLabel);
              oocxcoutD(static_cast<RooAbsArg*>(nullptr), InputArguments) << msgPrefix
                      << "assigning pdf " << compPdf->GetName() << "(member of " << citem.second.pdf->GetName()
                      << ") to super label " << superLabel << std::endl ;
            } else {
              oocoutW(nullptr, InputArguments) << msgPrefix << "WARNING: No p.d.f. associated with label "
                  << stype.second << " for component RooSimultaneous p.d.f " << citem.second.pdf->GetName()
                  << "which is associated with master index label " << citem.first << std::endl ;
            }
          }
        }
      }
    }
  }

  return out;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooSimultaneous::RooSimultaneous(const RooSimultaneous& other, const char* name) :
  RooAbsPdf(other,name),
  _plotCoefNormSet("!plotCoefNormSet",this,other._plotCoefNormSet),
  _plotCoefNormRange(other._plotCoefNormRange),
  _partIntMgr(other._partIntMgr,this),
  _indexCat("indexCat",this,other._indexCat),
  _numPdf(other._numPdf)
{
  // Copy proxy list
  for(auto* proxy : static_range_cast<RooRealProxy*>(other._pdfProxyList)) {
    _pdfProxyList.Add(new RooRealProxy(proxy->GetName(),this,*proxy)) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSimultaneous::~RooSimultaneous()
{
  _pdfProxyList.Delete() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the p.d.f associated with the given index category name

RooAbsPdf* RooSimultaneous::getPdf(RooStringView catName) const
{
  RooRealProxy* proxy = static_cast<RooRealProxy*>(_pdfProxyList.FindObject(catName.c_str()));
  return proxy ? static_cast<RooAbsPdf*>(proxy->absArg()) : nullptr;
}



////////////////////////////////////////////////////////////////////////////////
/// Associate given PDF with index category state label 'catLabel'.
/// The name state must be already defined in the index category.
///
/// RooSimultaneous can function without having a PDF associated
/// with every single state. The normalization in such cases is taken
/// from the number of registered PDFs, but getVal() will fail if
/// called for an unregistered index state.
///
/// PDFs may not overlap (i.e. share any variables) with the index category (function).
/// \param[in] pdf PDF to be added.
/// \param[in] catLabel Name of the category state to be associated to the PDF.
/// \return `true` in case of failure.

bool RooSimultaneous::addPdf(const RooAbsPdf& pdf, const char* catLabel)
{
  // PDFs cannot overlap with the index category
  if (pdf.dependsOn(_indexCat.arg())) {
    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName() << "): PDF '" << pdf.GetName()
           << "' overlaps with index category '" << _indexCat.arg().GetName() << "'."<< std::endl ;
    return true ;
  }

  // Each index state can only have one PDF associated with it
  if (_pdfProxyList.FindObject(catLabel)) {
    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName() << "): index state '"
           << catLabel << "' has already an associated PDF." << std::endl ;
    return true ;
  }

  const RooSimultaneous* simPdf = dynamic_cast<const RooSimultaneous*>(&pdf) ;
  if (simPdf) {

    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName()
           << ") ERROR: you cannot add a RooSimultaneous component to a RooSimultaneous using addPdf()."
           << " Use the constructor with RooArgList if input p.d.f.s or the map<string,RooAbsPdf&> instead." << std::endl ;
    return true ;

  } else {

    // Create a proxy named after the associated index state
    TObject* proxy = new RooRealProxy(catLabel,catLabel,this,const_cast<RooAbsPdf&>(pdf));
    _pdfProxyList.Add(proxy) ;
    _numPdf += 1 ;
  }

  return false ;
}

////////////////////////////////////////////////////////////////////////////////
/// Examine the pdf components and check if one of them can be extended or must be extended.
/// It is enough to have one component that can be extended or must be extended to return the flag in
/// the total simultaneous pdf.

RooAbsPdf::ExtendMode RooSimultaneous::extendMode() const
{
   bool anyCanExtend = false;

   for (auto *proxy : static_range_cast<RooRealProxy *>(_pdfProxyList)) {
      auto &pdf = static_cast<RooAbsPdf const&>(proxy->arg());
      if (pdf.mustBeExtended())
         return MustBeExtended;
      anyCanExtend |= pdf.canBeExtended();
   }
   return anyCanExtend ? CanBeExtended : CanNotBeExtended;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current value:
/// the value of the PDF associated with the current index category state

double RooSimultaneous::evaluate() const
{
   // Retrieve the proxy by index name
   RooRealProxy *proxy = static_cast<RooRealProxy *>(_pdfProxyList.FindObject(_indexCat.label()));
   if(!proxy) {
      return 0;
   }

   double nEvtTot = 1.0;
   double nEvtCat = 1.0;

   // Calculate relative weighting factor for sim-pdfs of all extendable components
   if (canBeExtended()) {

      nEvtTot = 0;
      nEvtCat = 0;

      for (auto *proxy2 : static_range_cast<RooRealProxy *>(_pdfProxyList)) {
         auto &pdf2 = static_cast<RooAbsPdf const &>(proxy2->arg());
         if(!pdf2.canBeExtended()) {
            // If one of the pdfs can't be expected, reset the normalization
            // factor to one and break out of the loop.
            nEvtTot = 1.0;
            nEvtCat = 1.0;
            break;
         }
         const double nEvt = pdf2.expectedEvents(_normSet);
         nEvtTot += nEvt;
         if (proxy == proxy2) {
            // Matching by proxy by pointer rather than pdfs, because it's
            // possible to have the same pdf used in different states.
            nEvtCat += nEvt;
         }
      }
   }
   double catFrac = nEvtCat / nEvtTot;

   // Return the selected PDF value, normalized by the relative number of
   // expected events if applicable.
   return *proxy * catFrac;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of expected events: If the index is in nset,
/// then return the sum of the expected events of all components,
/// otherwise return the number of expected events of the PDF
/// associated with the current index category state

double RooSimultaneous::expectedEvents(const RooArgSet* nset) const
{
  if (nset->contains(_indexCat.arg())) {

    double sum(0) ;

    for(auto * proxy : static_range_cast<RooRealProxy*>(_pdfProxyList)) {
      sum += (static_cast<RooAbsPdf*>(proxy->absArg()))->expectedEvents(nset) ;
    }

    return sum ;

  } else {

    // Retrieve the proxy by index name
    RooRealProxy* proxy = static_cast<RooRealProxy*>(_pdfProxyList.FindObject(_indexCat.label())) ;

    //assert(proxy!=0) ;
    if (proxy==nullptr) return 0 ;

    // Return the selected PDF value, normalized by the number of index states
    return (static_cast<RooAbsPdf*>(proxy->absArg()))->expectedEvents(nset);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward determination of analytical integration capabilities to component p.d.f.s
/// A unique code is assigned to the combined integration capabilities of all associated
/// p.d.f.s

Int_t RooSimultaneous::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars,
                      const RooArgSet* normSet, const char* rangeName) const
{
  // Declare that we can analytically integrate all requested observables
  analVars.add(allVars) ;

  // Retrieve (or create) the required partial integral list
  Int_t code ;

  // Check if this configuration was created before
  CacheElem* cache = static_cast<CacheElem*>(_partIntMgr.getObj(normSet,&analVars,nullptr,RooNameReg::ptr(rangeName))) ;
  if (cache) {
    code = _partIntMgr.lastIndex() ;
    return code+1 ;
  }
  cache = new CacheElem ;

  // Create the partial integral set for this request
  for(auto * proxy : static_range_cast<RooRealProxy*>(_pdfProxyList)) {
    cache->_partIntList.addOwned(std::unique_ptr<RooAbsReal>{proxy->arg().createIntegral(analVars,normSet,nullptr,rangeName)});
  }

  // Store the partial integral list and return the assigned code ;
  code = _partIntMgr.setObj(normSet,&analVars,cache,RooNameReg::ptr(rangeName)) ;

  return code+1 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return analytical integration defined by given code

double RooSimultaneous::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* /*rangeName*/) const
{
  // No integration scenario
  if (code==0) {
    return getVal(normSet) ;
  }

  // Partial integration scenarios, rangeName already encoded in 'code'
  CacheElem* cache = static_cast<CacheElem*>(_partIntMgr.getObjByIndex(code-1)) ;

  RooRealProxy* proxy = static_cast<RooRealProxy*>(_pdfProxyList.FindObject(_indexCat.label())) ;
  Int_t idx = _pdfProxyList.IndexOf(proxy) ;
  return (static_cast<RooAbsReal*>(cache->_partIntList.at(idx)))->getVal(normSet) ;
}






////////////////////////////////////////////////////////////////////////////////
/// Back-end for plotOn() implementation on RooSimultaneous which
/// needs special handling because a RooSimultaneous PDF cannot
/// project out its index category via integration. plotOn() will
/// abort if this is requested without providing a projection dataset.

RooPlot* RooSimultaneous::plotOn(RooPlot *frame, RooLinkedList& cmdList) const
{
  // Sanity checks
  if (plotSanityChecks(frame)) return frame ;

  // Special case: if an asymmetry is requested with respect to our index
  // category, we cannot reroute the plotting to the component pdfs. The
  // component pdfs don't depend on the index category, so the asymmetry engine
  // in the base class would not be able to split them by index state. Instead,
  // we delegate directly to the base class implementation, which constructs the
  // asymmetry from the two index-state component pdfs (see the overridden
  // createAsymmetryComponent() and GitHub issue #14255).
  if (auto *asymCmd = static_cast<RooCmdArg *>(cmdList.FindObject("Asymmetry"))) {
    auto *asymCat = dynamic_cast<RooAbsCategory const *>(asymCmd->getObject(0));
    if (asymCat && asymCat == &_indexCat.arg()) {

      RooLinkedList cmdList2(cmdList);

      // The base-class asymmetry-plotting engine averages the projection over
      // the projection dataset. This is not supported for the composite data
      // stores that back datasets with a category index, so we flatten such a
      // projection dataset into a plain (vector-backed) copy first. Both the
      // copy and the replacement command must outlive the plotOn() call below,
      // because the command list only stores pointers to them.
      std::unique_ptr<RooAbsData> flatProjData;
      RooCmdArg newProjWData;
      if (auto *projWData = static_cast<RooCmdArg *>(cmdList2.FindObject("ProjData"))) {
        auto *projData = dynamic_cast<RooDataSet const *>(projWData->getObject(1));
        if (projData && dynamic_cast<RooCompositeDataStore const *>(projData->store())) {
          flatProjData = std::make_unique<RooDataSet>(projData->GetName(), projData->GetTitle(), *projData->get(),
                                                      RooFit::Import(*const_cast<RooDataSet *>(projData)));
          const RooArgSet *projDataSet = projWData->getSet(0);
          newProjWData = projDataSet ? RooFit::ProjWData(*projDataSet, *flatProjData)
                                     : RooFit::ProjWData(*flatProjData);
          replaceOrAdd(cmdList2, newProjWData);
        }
      }

      return RooAbsReal::plotOn(frame, cmdList2);
    }
  }

  // Extract projection configuration from command list
  RooCmdConfig pc("RooSimultaneous::plotOn(" + std::string(GetName()) + ")");
  pc.defineString("sliceCatState","SliceCat",0,"",true) ;
  pc.defineDouble("scaleFactor","Normalization",0,1.0) ;
  pc.defineInt("scaleType","Normalization",0,RooAbsPdf::Relative) ;
  pc.defineObject("sliceCatList","SliceCat",0,nullptr,true) ;
  // This dummy is needed for plotOn to recognize the "SliceCatMany" command.
  // It is not used directly, but the "SliceCat" commands are nested in it.
  // Removing this dummy definition results in "ERROR: unrecognized command: SliceCatMany".
  pc.defineObject("dummy1","SliceCatMany",0) ;
  pc.defineSet("projSet","Project",0) ;
  pc.defineSet("sliceSet","SliceVars",0) ;
  pc.defineSet("projDataSet","ProjData",0) ;
  pc.defineObject("projData","ProjData",1) ;
  pc.defineMutex("Project","SliceVars") ;
  pc.allowUndefined() ; // there may be commands we don't handle here

  // Process and check varargs
  pc.process(cmdList) ;
  if (!pc.ok(true)) {
    return frame ;
  }

  RooAbsData* projData = static_cast<RooAbsData*>(pc.getObject("projData")) ;
  const RooArgSet* projDataSet = pc.getSet("projDataSet");
  const RooArgSet* sliceSetTmp = pc.getSet("sliceSet") ;
  std::unique_ptr<RooArgSet> sliceSet( sliceSetTmp ? (static_cast<RooArgSet*>(sliceSetTmp->Clone())) : nullptr );
  const RooArgSet* projSet = pc.getSet("projSet") ;
  double scaleFactor = pc.getDouble("scaleFactor") ;
  ScaleType stype = (ScaleType) pc.getInt("scaleType") ;


  // Look for category slice arguments and add them to the master slice list if found
  const char* sliceCatState = pc.getString("sliceCatState",nullptr,true) ;
  const RooLinkedList& sliceCatList = pc.getObjectList("sliceCatList") ;
  if (sliceCatState) {

    // Make the master slice set if it doesnt exist
    if (!sliceSet) {
      sliceSet = std::make_unique<RooArgSet>();
    }

    // Prepare comma separated label list for parsing
    auto catTokens = ROOT::Split(sliceCatState, ",");

    // Loop over all categories provided by (multiple) Slice() arguments
    unsigned int tokenIndex = 0;
    for(auto * scat : static_range_cast<RooCategory*>(sliceCatList)) {
      const char* slabel = tokenIndex >= catTokens.size() ? nullptr : catTokens[tokenIndex++].c_str();

      if (slabel) {
        // Set the slice position to the value indicated by slabel
        scat->setLabel(slabel) ;
        // Add the slice category to the master slice set
        sliceSet->add(*scat,false) ;
      }
    }
  }

  // Check if we have a projection dataset
  if (!projData) {
    coutE(InputArguments) << "RooSimultaneous::plotOn(" << GetName() << ") ERROR: must have a projection dataset for index category" << std::endl ;
    return frame ;
  }

  // Make list of variables to be projected
  RooArgSet projectedVars ;
  if (sliceSet) {
    makeProjectionSet(frame->getPlotVar(),frame->getNormVars(),projectedVars,true) ;

    // Take out the sliced variables
    for (const auto sliceArg : *sliceSet) {
      RooAbsArg* arg = projectedVars.find(sliceArg->GetName()) ;
      if (arg) {
        projectedVars.remove(*arg) ;
      } else {
        coutI(Plotting) << "RooAbsReal::plotOn(" << GetName() << ") slice variable "
            << sliceArg->GetName() << " was not projected anyway" << std::endl ;
      }
    }
  } else if (projSet) {
    makeProjectionSet(frame->getPlotVar(),projSet,projectedVars,false) ;
  } else {
    makeProjectionSet(frame->getPlotVar(),frame->getNormVars(),projectedVars,true) ;
  }

  bool projIndex(false) ;

  if (!_indexCat.arg().isDerived()) {
    // *** Error checking for a fundamental index category ***
    //cout << "RooSim::plotOn: index is fundamental" << std::endl ;

    // Check that the provided projection dataset contains our index variable
    if (!projData->get()->find(_indexCat.arg().GetName())) {
      coutE(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") ERROR: Projection over index category "
            << "requested, but projection data set doesn't contain index category" << std::endl ;
      return frame ;
    }

    if (projectedVars.find(_indexCat.arg().GetName())) {
      projIndex=true ;
    }

  } else {
    // *** Error checking for a composite index category ***

    // Determine if any servers of the index category are in the projectedVars
    RooArgSet projIdxServers ;
    bool anyServers(false) ;
    for (const auto server : flattenedCatList()) {
      if (projectedVars.find(server->GetName())) {
        anyServers=true ;
        projIdxServers.add(*server) ;
      }
    }

    // Check that the projection dataset contains all the
    // index category components we're projecting over

    // Determine if all projected servers of the index category are in the projection dataset
    bool allServers(true) ;
    std::string missing;
    for (const auto server : projIdxServers) {
      if (!projData->get()->find(server->GetName())) {
        allServers=false ;
        missing = server->GetName();
      }
    }

    if (!allServers) {
      coutE(Plotting) << "RooSimultaneous::plotOn(" << GetName()
          << ") ERROR: Projection dataset doesn't contain complete set of index categories to do projection."
          << "\n\tcategory " << missing << " is missing." << std::endl ;
      return frame ;
    }

    if (anyServers) {
      projIndex = true ;
    }
  }

  // Calculate relative weight fractions of components
  std::unique_ptr<Roo1DTable> wTable( projData->table(_indexCat.arg()) );

  // Clone the index category to be able to cycle through the category states for plotting without
  // affecting the category state of our instance
  RooArgSet idxCloneSet;
  RooArgSet(*_indexCat).snapshot(idxCloneSet, true);
  auto idxCatClone = static_cast<RooAbsCategoryLValue*>(idxCloneSet.find(_indexCat->GetName()) );
  assert(idxCatClone);

  // Make list of category columns to exclude from projection data
  std::unique_ptr<RooArgSet> idxCompSliceSet( idxCatClone->getObservables(frame->getNormVars()) );

  // If we don't project over the index, just do the regular plotOn
  if (!projIndex) {

    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
          << " represents a slice in the index category ("  << _indexCat.arg().GetName() << ")" << std::endl ;

    // Reduce projData: take out fitCat (component) columns and entries that don't match selected slice
    // Construct cut string to only select projection data event that match the current slice

    // Make cut string to exclude rows from projection data
    if (sliceSet) {
      for (auto *idxComp : static_range_cast<RooCategory *>(*idxCompSliceSet)) {
        if (auto* slicedComponent = static_cast<const RooAbsCategory*>(sliceSet->find(*idxComp))) {
          idxComp->setIndex(slicedComponent->getCurrentIndex(), false);
        }
      }
    }
    std::string cutString = RooFit::Detail::makeSliceCutString(*idxCompSliceSet);

    // Make temporary projData without RooSim index category components
    RooArgSet projDataVars(*projData->get()) ;
    projDataVars.remove(*idxCompSliceSet,true,true) ;

    std::unique_ptr<RooAbsData>
       projDataTmp(projData->reduce(RooFit::SelectVars(projDataVars), RooFit::Cut(cutString.c_str())));

    // Override normalization and projection dataset
    RooCmdArg tmp1 =
       RooFit::Normalization(scaleFactor * wTable->get(idxCatClone->getCurrentLabel()), RooAbsReal::NumEvent);
    RooCmdArg tmp2 = RooFit::ProjWData(*projDataSet,*projDataTmp) ;

    // WVE -- do not adjust normalization for asymmetry plots
    RooLinkedList cmdList2(cmdList) ;
    if (!cmdList.find("Asymmetry")) {
      replaceOrAdd(cmdList2, tmp1);
    }
    replaceOrAdd(cmdList2, tmp2);

    // Plot single component
    RooPlot* retFrame = getPdf(idxCatClone->getCurrentLabel())->plotOn(frame,cmdList2);
    return retFrame ;
  }

  // If we project over the index, plot using a temporary RooAddPdf
  // using the weights from the data as coefficients

  // Build the list of indexCat components that are sliced
  idxCompSliceSet->remove(projectedVars,true,true) ;

  // Make a new expression that is the weighted sum of requested components
  RooArgList pdfCompList ;
  RooArgList wgtCompList ;
//RooAbsPdf* pdf ;
  double sumWeight(0) ;
  for(auto * proxy : static_range_cast<RooRealProxy*>(_pdfProxyList)) {

    idxCatClone->setLabel(proxy->name()) ;

    // Determine if this component is the current slice (if we slice)
    bool skip(false) ;
    for (const auto idxSliceCompArg : *idxCompSliceSet) {
      const auto idxSliceComp = static_cast<RooAbsCategory*>(idxSliceCompArg);
      RooAbsCategory* idxComp = static_cast<RooAbsCategory*>(idxCloneSet.find(idxSliceComp->GetName())) ;
      if (idxComp->getCurrentIndex()!=idxSliceComp->getCurrentIndex()) {
        skip=true ;
        break ;
      }
    }
    if (skip) continue ;

    // Instantiate a RRV holding this pdfs weight
    wgtCompList.addOwned(std::make_unique<RooRealVar>(proxy->name(),"coef",wTable->get(proxy->name())));
    sumWeight += wTable->getFrac(proxy->name()) ;

    // Add the PDF to list list
    pdfCompList.add(proxy->arg()) ;
  }

  TString plotVarName(GetName()) ;
  RooAddPdf plotVar{plotVarName,"weighted sum of RS components",pdfCompList,wgtCompList};

  // Fix appropriate coefficient normalization in plot function
  if (!_plotCoefNormSet.empty()) {
    plotVar.fixAddCoefNormalization(_plotCoefNormSet) ;
  }

  std::unique_ptr<RooAbsData> projDataTmp;
  RooArgSet projSetTmp ;
  if (projData) {

    // Construct cut string to only select projection data event that match the current slice
    std::string cutString = RooFit::Detail::makeSliceCutString(*idxCompSliceSet);

    // Make temporary projData without RooSim index category components
    RooArgSet projDataVars(*projData->get()) ;
    RooArgSet idxCatServers;
    _indexCat.arg().getObservables(frame->getNormVars(), idxCatServers) ;

    projDataVars.remove(idxCatServers,true,true) ;

    projDataTmp = std::unique_ptr<RooAbsData>{projData->reduce(RooFit::SelectVars(projDataVars), RooFit::Cut(cutString.c_str()))};



    if (projSet) {
      projSetTmp.add(*projSet) ;
      projSetTmp.remove(idxCatServers,true,true);
    }
  }


  if (_indexCat.arg().isDerived() && !idxCompSliceSet->empty()) {
    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
          << " represents a slice in index category components " << *idxCompSliceSet << std::endl ;

    RooArgSet idxCompProjSet;
    _indexCat.arg().getObservables(frame->getNormVars(), idxCompProjSet) ;
    idxCompProjSet.remove(*idxCompSliceSet,true,true) ;
    if (!idxCompProjSet.empty()) {
      coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
            << " averages with data index category components " << idxCompProjSet << std::endl ;
    }
  } else {
    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
          << " averages with data index category (" << _indexCat.arg().GetName() << ")" << std::endl ;
  }


  // Override normalization and projection dataset
  RooLinkedList cmdList2(cmdList) ;

  RooCmdArg tmp1 = RooFit::Normalization(scaleFactor*sumWeight,stype) ;
  RooCmdArg tmp2 = RooFit::ProjWData(*projDataSet,*projDataTmp) ;
  // WVE -- do not adjust normalization for asymmetry plots
  if (!cmdList.find("Asymmetry")) {
    replaceOrAdd(cmdList2, tmp1);
  }
  replaceOrAdd(cmdList2, tmp2);

  RooPlot* frame2 ;
  if (!projSetTmp.empty()) {
    // Plot temporary function
    RooCmdArg tmp3 = RooFit::Project(projSetTmp) ;
    replaceOrAdd(cmdList2, tmp3);
    frame2 = plotVar.plotOn(frame,cmdList2) ;
  } else {
    // Plot temporary function
    frame2 = plotVar.plotOn(frame,cmdList2) ;
  }

  return frame2 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Build the component function of an asymmetry plot (see
/// RooAbsReal::plotAsymOn()) for a fixed state of the asymmetry category.
///
/// When the asymmetry is requested in our own index category, the component for
/// a given index state is simply the corresponding pdf. We return a clone of
/// that pdf directly instead of a RooSimultaneous with a pinned index, because
/// a RooSimultaneous compiles its per-category observables with a category
/// prefix. That prefix makes it incompatible with the vectorized evaluation
/// backend that averages the asymmetry over the projection data, and would
/// otherwise silently yield a flat (zero) asymmetry (see issue #14255). For any
/// other asymmetry category we fall back to the generic implementation.

std::unique_ptr<RooAbsReal>
RooSimultaneous::createAsymmetryComponent(const RooAbsCategoryLValue &asymCat, const RooAbsCategoryLValue &asymCatState) const
{
   if (&asymCat == &_indexCat.arg()) {
      const std::string &label = _indexCat.arg().lookupName(asymCatState.getCurrentIndex());
      if (RooAbsPdf *pdf = getPdf(label)) {
         return RooHelpers::cloneTreeWithSameParameters(static_cast<RooAbsReal const &>(*pdf));
      }
   }
   return RooAbsReal::createAsymmetryComponent(asymCat, asymCatState);
}


////////////////////////////////////////////////////////////////////////////////
/// Interface function used by test statistics to freeze choice of observables
/// for interpretation of fraction coefficients. Needed here because a RooSimultaneous
/// works like a RooAddPdf when plotted

void RooSimultaneous::selectNormalization(const RooArgSet* normSet, bool /*force*/)
{
  _plotCoefNormSet.removeAll() ;
  if (normSet) {
     // The index category must not be stored in the set: it is meaningless for
     // the coefficient normalization, since it is never an observable of the
     // component pdfs (RooAddPdf::selectNormalization() would filter it out
     // again anyway). Worse, it is already registered as a value server via
     // the index category proxy, and registering the same server a second
     // time through this non-propagating set proxy corrupts the reference
     // counts of the server's client lists when the set is cleared again.
     RooArgSet filteredNormSet{*normSet};
     filteredNormSet.remove(_indexCat.arg(), true, true);
     _plotCoefNormSet.add(filteredNormSet);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Interface function used by test statistics to freeze choice of range
/// for interpretation of fraction coefficients. Needed here because a RooSimultaneous
/// works like a RooAddPdf when plotted

void RooSimultaneous::selectNormalizationRange(const char* normRange2, bool /*force*/)
{
  _plotCoefNormRange = RooNameReg::ptr(normRange2) ;
}




////////////////////////////////////////////////////////////////////////////////

RooAbsGenContext* RooSimultaneous::autoGenContext(const RooArgSet &vars, const RooDataSet* prototype,
                    const RooArgSet* auxProto, bool verbose, bool autoBinned, const char* binnedTag) const
{
  const char* idxCatName = _indexCat.arg().GetName() ;

  if (vars.find(idxCatName) && prototype==nullptr
      && (auxProto==nullptr || auxProto->empty())
      && (autoBinned || (binnedTag && strlen(binnedTag)))) {

    // Return special generator config that can also do binned generation for selected states
    return new RooSimSplitGenContext(*this,vars,verbose,autoBinned,binnedTag) ;

  } else {

    // Return regular generator config ;
    return genContext(vars,prototype,auxProto,verbose) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Return specialized generator context for simultaneous p.d.f.s

RooAbsGenContext* RooSimultaneous::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                     const RooArgSet* auxProto, bool verbose) const
{
  RooArgSet allVars{vars};
  if(prototype) allVars.add(*prototype->get());

  RooArgSet catsAmongAllVars;
  allVars.selectCommon(flattenedCatList(), catsAmongAllVars);

  // Not generating index cat: we better error out because it's not clear what
  // the user expects here. Does the user want to generate according to the
  // currently-selected pdf? Or does the user want to generate global
  // observable values according to the union of all category pdfs?
  // Print an error and tell the user what to do to explicitly.
  if(catsAmongAllVars.empty()) {
    coutE(InputArguments) << "RooSimultaneous::generateSimGlobal(" << GetName()
           << ") asking to generate without the index category!\n"
           << "It's not clear what to do. you probably want to either:\n"
           << "\n"
           << "    1. Generate according to the currently-selected pdf.\n"
           << "       Please do this explicitly with:\n"
           << "           simpdf->getPdf(simpdf->indexCat().getCurrentLabel())->generate(vars, ...)\n"
           << "\n"
           << "    1. Generate global observable values according to the union of all component pdfs.\n"
           << "       For this, please use simpdf->generateSimGlobal(vars, ...)\n"
           << std::endl;
    return nullptr;
  }

  RooArgSet catsAmongProtoVars;
  if(prototype) {
    prototype->get()->selectCommon(flattenedCatList(), catsAmongProtoVars);

    if(!catsAmongProtoVars.empty() && catsAmongProtoVars.size() != flattenedCatList().size()) {
      // Abort if we have only part of the servers
      coutE(Plotting) << "RooSimultaneous::genContext: ERROR: prototype must include either all "
            << " components of the RooSimultaneous index category or none " << std::endl;
      return nullptr;
    }
  }

  return new RooSimGenContext(*this,vars,prototype,auxProto,verbose) ;
}




////////////////////////////////////////////////////////////////////////////////

RooDataHist* RooSimultaneous::fillDataHist(RooDataHist *hist,
                                           const RooArgSet* nset,
                                           double scaleFactor,
                                           bool correctForBinVolume,
                                           bool showProgress) const
{
  if (RooAbsReal::fillDataHist (hist, nset, scaleFactor,
                                correctForBinVolume, showProgress) == nullptr)
    return nullptr;

  const double sum = hist->sumEntries();
  if (sum != 0) {
    for (int i=0 ; i<hist->numEntries() ; i++) {
      hist->set(i, hist->weight(i) / sum, 0.);
    }
  }

  return hist;
}




////////////////////////////////////////////////////////////////////////////////
/// Special generator interface for generation of 'global observables' -- for RooStats tools.
///
/// \note Why one can't just use RooAbsPdf::generate()? That's because when
/// using the regular generate() method, a specific component pdf is selected
/// for each generated entry according to the index category value. However,
/// global observable values are independent of the current index category,
/// which can best be illustrated with the case where a global observable
/// corresponds to a nuisance parameter that is relevant for multiple channels.
/// So the interpretation of what is an entry in the generated dataset is very
/// different, hence the separate function.

RooFit::OwningPtr<RooDataSet> RooSimultaneous::generateSimGlobal(const RooArgSet& whatVars, Int_t nEvents)
{
  // Generating the index category together with the global observables doesn't make any sense.
  RooArgSet catsAmongAllVars;
  whatVars.selectCommon(flattenedCatList(), catsAmongAllVars);
  if(!catsAmongAllVars.empty()) {
    coutE(InputArguments) << "RooSimultaneous::generateSimGlobal(" << GetName()
           << ") asking to generate global obserables at the same time as the index category!\n"
           << "This doesn't make any sense: global observables are generally not related to a specific channel.\n"
           << std::endl;
    return nullptr;
  }

  // Make set with clone of variables (placeholder for output)
  RooArgSet globClone;
  whatVars.snapshot(globClone);

  auto data = std::make_unique<RooDataSet>("gensimglobal","gensimglobal",whatVars);

  for (Int_t i=0 ; i<nEvents ; i++) {
    for (const auto& nameIdx : indexCat()) {

      // Get pdf associated with state from simpdf
      RooAbsPdf* pdftmp = getPdf(nameIdx.first);

      RooArgSet globtmp;
      pdftmp->getObservables(&whatVars, globtmp) ;

      // If there are any, generate only global variables defined by the pdf
      // associated with this state and transfer values to output placeholder.
      if (!globtmp.empty()) {
        globClone.assign(*std::unique_ptr<RooDataSet>{pdftmp->generate(globtmp,1)}->get(0)) ;
      }
    }
    data->add(globClone) ;
  }

  return RooFit::makeOwningPtr(std::move(data));
}


/// Wraps the components of this RooSimultaneous in RooBinSamplingPdfs.
/// \param[in] data The dataset to be used in the eventual fit, used to figure
///            out the observables and whether the dataset is binned.
/// \param[in] precision Precision argument for all created RooBinSamplingPdfs.
void RooSimultaneous::wrapPdfsInBinSamplingPdfs(RooAbsData const &data, double precision) {

  if (precision < 0.) return;

  RooArgSet newSamplingPdfs;

  for (auto const &item : this->indexCat()) {

    auto const &catName = item.first;
    auto &pdf = *this->getPdf(catName);

    if (auto newSamplingPdf = RooBinSamplingPdf::create(pdf, data, precision)) {
      // Set the "ORIGNAME" attribute the indicate to
      // RooAbsArg::redirectServers() which pdf should be replaced by this
      // RooBinSamplingPdf in the RooSimultaneous.
      newSamplingPdf->setAttribute(
          (std::string("ORIGNAME:") + pdf.GetName()).c_str());
      newSamplingPdfs.addOwned(std::move(newSamplingPdf));
    }
  }

  this->redirectServers(newSamplingPdfs, false, true);
  this->addOwnedComponents(std::move(newSamplingPdfs));
}


/// Wraps the components of this RooSimultaneous in RooBinSamplingPdfs, with a
/// different precision parameter for each component.
/// \param[in] data The dataset to be used in the eventual fit, used to figure
///            out the observables and whether the dataset is binned.
/// \param[in] precisions The map that gives the precision argument for each
///            component in the RooSimultaneous. The keys are the pdf names. If
///            there is no value for a given component, it will not use the bin
///            integration. Otherwise, the value has the same meaning than in
///            the IntegrateBins() command argument for RooAbsPdf::fitTo().
/// \param[in] useCategoryNames If this flag is set, the category names will be
///            used to look up the precision in the precisions map instead of
///            the pdf names.
void RooSimultaneous::wrapPdfsInBinSamplingPdfs(RooAbsData const &data,
                                                std::map<std::string, double> const& precisions,
                                                bool useCategoryNames /*=false*/) {

  constexpr double defaultPrecision = -1.;

  RooArgSet newSamplingPdfs;

  for (auto const &item : this->indexCat()) {

    auto const &catName = item.first;
    auto &pdf = *this->getPdf(catName);
    std::string pdfName = pdf.GetName();

    auto found = precisions.find(useCategoryNames ? catName : pdfName);
    const double precision =
        found != precisions.end() ? found->second : defaultPrecision;
    if (precision < 0.)
      continue;

    if (auto newSamplingPdf = RooBinSamplingPdf::create(pdf, data, precision)) {
      // Set the "ORIGNAME" attribute the indicate to
      // RooAbsArg::redirectServers() which pdf should be replaced by this
      // RooBinSamplingPdf in the RooSimultaneous.
      newSamplingPdf->setAttribute(
          (std::string("ORIGNAME:") + pdf.GetName()).c_str());
      newSamplingPdfs.addOwned(std::move(newSamplingPdf));
    }
  }

  this->redirectServers(newSamplingPdfs, false, true);
  this->addOwnedComponents(std::move(newSamplingPdfs));
}

/// Internal utility function to get a list of all category components for this
/// RooSimultaneous. The output contains only the index category if it is a
/// RooCategory, or the list of all category components if it is a
/// RooSuperCategory.
RooArgSet const& RooSimultaneous::flattenedCatList() const
{
   // Note that the index category of a RooSimultaneous can only be of type
   // RooCategory or RooSuperCategory, because these are the only classes that
   // inherit from RooAbsCategoryLValue.
   if (auto superCat = dynamic_cast<RooSuperCategory const*>(&_indexCat.arg())) {
       return superCat->inputCatList();
   }

   if(!_indexCatSet) {
      _indexCatSet = std::make_unique<RooArgSet>(_indexCat.arg());
   }
   return *_indexCatSet;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the index category is among the variables `vars`, matching by
/// name. For a RooSuperCategory index, any of its input categories counts.
///
/// If the index category is not among the observables of a fit, this
/// RooSimultaneous does not split the data into channels: it acts as a
/// "switch" that evaluates to the component selected by the current index
/// state, analogous to RooMultiPdf. Fitting infrastructure uses this check to
/// decide between the two modes.
bool RooSimultaneous::indexCatIsObservable(RooArgSet const &vars) const
{
   RooArgSet catsAmongVars;
   vars.selectCommon(flattenedCatList(), catsAmongVars);
   return !catsAmongVars.empty();
}

namespace {

/// Whether the experimental compilation of a RooSimultaneous into an ordinary
/// mixture pdf is requested via the ROOFIT_SIM_COMPILE_MIXTURE environment
/// variable.
bool simMixtureCompileRequested()
{
   const char *env = std::getenv("ROOFIT_SIM_COMPILE_MIXTURE");
   return env && *env && std::string_view{env} != "0";
}

/// Variant of the mixture compilation for a simultaneous pdf whose channels
/// all use the binned likelihood optimization. The compiled pdf is a single
/// unnormalized sum of the indicator-gated channel yields,
/// \f[
///   Y(\vec{x}, c) = \sum_s \mathbf{1}[c = s] \; Y_s(\vec{x}),
/// \f]
/// built as a RooRealSumPdf of RooProducts with unit coefficients. It is
/// compiled through the standard binned-likelihood machinery of
/// RooRealSumPdf::compileForNormSet(): the RooBinWidthFunctions in the
/// channel pdfs disable themselves, so the compiled values are directly the
/// expected bin yields, and the compiled pdf carries the
/// "BinnedLikelihoodActive(Yields)" attributes that make RooNLLVarNew sum
/// Poisson terms over the concatenated bins and make the data loading retain
/// zero-weight entries. The result is identical to the sum of the per-channel
/// binned likelihoods of the channel-splitting path; the "SimCount" attribute
/// reproduces the legacy convention of adding sumOfWeights * log(nChannels).
template <typename FallBackFunc>
std::unique_ptr<RooAbsArg>
compileSimPdfAsBinnedMixture(RooSimultaneous const &simPdf, RooArgSet const &normSet,
                             RooFit::Detail::CompileContext &ctx, FallBackFunc const &fallBack)
{
   RooAbsCategoryLValue const &indexCat = simPdf.indexCat();

   // The stand-in for the index category. The range is set once all state
   // indices are known.
   auto standIn = std::make_unique<RooRealVar>(indexCat.GetName(), indexCat.GetTitle(), 0.0);

   RooArgList prods;
   int minIndex = std::numeric_limits<int>::max();
   int maxIndex = std::numeric_limits<int>::min();

   std::map<std::string, RooAbsPdf const *> seenChannelPdfs;

   for (auto const &catState : indexCat) {
      RooAbsPdf *channelPdf = simPdf.getPdf(catState.first.c_str());

      auto [seenIt, inserted] = seenChannelPdfs.emplace(channelPdf->GetName(), channelPdf);
      if (!inserted && seenIt->second != channelPdf) {
         // Two different pdf objects with the same name would collide in the
         // name-keyed deduplication of the graph compilation. Attaching the
         // same pdf object to several channels is fine here, unlike in the
         // extended unbinned case: the binned likelihood involves no
         // per-channel expected-events functions that could collide.
         return fallBack("two channels use different pdfs with the same name \"" + seenIt->first + "\"");
      }

      // The concatenated binned likelihood can only interpret the compiled
      // pdf values directly as yields (the "BinnedLikelihoodActiveYields"
      // mode). That requires RooBinWidthFunctions in the channel pdfs, which
      // disable themselves during the binned-likelihood compilation. Without
      // them, the likelihood would have to multiply by per-channel bin
      // volumes, which the single concatenated RooNLLVarNew doesn't support.
      RooAbsPdf const &binnedPdf = *RooHelpers::getBinnedL(*channelPdf).binnedPdf;
      RooArgList binnedPdfNodes;
      binnedPdf.treeNodeServerList(&binnedPdfNodes);
      bool hasBinWidthFunc = false;
      for (RooAbsArg const *node : binnedPdfNodes) {
         if (dynamic_cast<RooBinWidthFunction const *>(node)) {
            hasBinWidthFunc = true;
            break;
         }
      }
      if (!hasBinWidthFunc) {
         return fallBack("the binned-likelihood pdf of channel \"" + catState.first +
                         "\" has no RooBinWidthFunction, so its values can't be interpreted as bin yields");
      }

      minIndex = std::min(minIndex, catState.second);
      maxIndex = std::max(maxIndex, catState.second);

      // The suffixes make the names collision-safe, see the unbinned mixture
      // compilation below.
      std::string baseName = std::string(simPdf.GetName()) + "_" + catState.first;
      auto indicator = std::make_unique<RooFit::Detail::RooChannelIndicatorPdf>(
         (baseName + "_mixtureIndicator").c_str(), (baseName + "_mixtureIndicator").c_str(), *standIn, catState.second);
      // Declare to the RooFit::Evaluator that this node is a data-only
      // {0,1}-valued mask, so that it can restrict the evaluation of the
      // other factors in the gated product to the events of this channel (see
      // Evaluator::rangeRestrictionAnalysis()). Without this, every channel
      // pdf is evaluated on the concatenated bins of all channels.
      indicator->setAttribute("BinaryMask");
      auto prod = std::make_unique<RooProduct>((baseName + "_mixtureTerm").c_str(), (baseName + "_mixtureTerm").c_str(),
                                               RooArgList(*indicator, *channelPdf));
      prod->addOwnedComponents(std::move(indicator));
      prods.addOwned(std::move(prod));
   }

   const std::size_t nChannels = prods.size();

   standIn->setRange(minIndex - 0.5, maxIndex + 0.5);
   standIn->setVal(indexCat.getCurrentIndex());

   // Unit coefficients for the sum of the gated channel yields. An owned
   // constant is used instead of RooFit::RooConst(), because the global
   // constants registry must not end up in a compiled computation graph
   // (concurrent evaluators would clash on its data token).
   std::string coefName = std::string(simPdf.GetName()) + "_mixtureCoef";
   auto coefVar = std::make_unique<RooConstVar>(coefName.c_str(), coefName.c_str(), 1.0);
   RooArgList coefs;
   for (std::size_t i = 0; i < nChannels; ++i) {
      coefs.add(*coefVar);
   }

   auto mixture = std::make_unique<RooRealSumPdf>(simPdf.GetName(), simPdf.GetTitle(), prods, coefs);
   // Request the binned-likelihood compilation of RooRealSumPdf for the
   // mixture sum itself.
   mixture->setAttribute("BinnedLikelihood");

   // The normalization set for the mixture, with the index category replaced
   // by its real-valued stand-in.
   RooArgSet mixtureNormSet;
   for (RooAbsArg *arg : normSet) {
      mixtureNormSet.add(arg->namePtr() == indexCat.namePtr() ? *standIn : *arg);
   }

   mixture->addOwnedComponents(std::move(standIn));
   mixture->addOwnedComponents(std::move(prods));
   mixture->addOwnedComponents(std::move(coefVar));

   std::unique_ptr<RooAbsArg> compiled = mixture->compileForNormSet(mixtureNormSet, ctx);

   if (!compiled->getAttribute("BinnedLikelihoodActiveYields")) {
      // The RooBinWidthFunction check above should have guaranteed the yields
      // mode; without it, the Poisson terms would silently use probability
      // densities as yields.
      throw std::runtime_error("RooSimultaneous::compileForNormSet(): the binned-likelihood mixture compilation "
                               "unexpectedly didn't end up in yields mode");
   }

   // The per-channel binned likelihoods of the channel-splitting path each
   // add the legacy sumOfWeights * log(nChannels) term, so the concatenated
   // likelihood has to be asked to do the same.
   compiled->setStringAttribute("SimCount", std::to_string(nChannels).c_str());

   // Mark the compiled mixture terms as products gated by a binary mask, so
   // that the RooFit::Evaluator can restrict the evaluation of the channel
   // pdfs to the bins of their own channel. The compiled product nodes are
   // new objects that don't inherit the attributes of the RooProducts above,
   // so the marking has to happen after the compilation. The sum over the
   // gated products is exact for the skipped bins too: their buffer entries
   // are exact zeros.
   if (auto *compiledSumPdf = dynamic_cast<RooRealSumPdf *>(compiled.get())) {
      for (RooAbsArg *component : compiledSumPdf->funcList()) {
         for (RooAbsArg *server : component->servers()) {
            if (server->getAttribute("BinaryMask") && server->isValueServer(*component)) {
               component->setAttribute("MaskGatedProduct");
               break;
            }
         }
      }
   }

   // Keep the uncompiled mixture template alive: some normalization sets
   // stored inside RooProdPdf are disconnected from the computation graph, so
   // server redirection has no control over them (see the comment in the
   // channel-splitting compilation below). Rename it to avoid a name clash
   // with the compiled pdf.
   mixture->SetName((std::string("_") + mixture->GetName()).c_str());
   compiled->addOwnedComponents(std::move(mixture));

   return compiled;
}

/// Experimental alternative to the channel-splitting compilation below: the
/// RooSimultaneous replaces itself with an ordinary mixture pdf,
/// \f[
///   P(\vec{x}, c) = \sum_s w_s \; \mathbf{1}[c = s] \; \mathrm{pdf}_s(\vec{x}),
/// \f]
/// built from standard components (RooAddPdf, RooProdPdf and the indicator
/// densities), so that everything downstream in the likelihood creation --
/// RooNLLVarNew, RooEvaluatorWrapper, RooFit::Evaluator and the data loading
/// -- can stay completely agnostic of the simultaneous structure. The index
/// category is replaced by a real-valued stand-in variable with the same
/// name, which gets filled directly from the (double-converted) category
/// column of the dataset. No category node is left in the compiled graph.
///
/// The mixture weights are constant \f$ w_s = 1/C \f$ for a non-extended fit,
/// which reproduces the legacy convention of adding
/// \f$ \sum_i w_i \log(C) \f$ to the NLL (with \f$ C \f$ the number of
/// channels) exactly. For an extended fit, the RooAddPdf is built in
/// all-extendable mode, where the coefficients are the expected event yields
/// of the channels: the resulting single extended NLL is then equal to the
/// sum of the per-channel extended NLLs.
///
/// Note that no explicit "padding" pdfs over the observables that a channel
/// does not depend on are needed: the factorized normalization of RooProdPdf
/// normalizes each factor over its own observables only, which is
/// mathematically identical to padding each channel with uniform densities
/// over the unused observables and dividing out their constant volumes.
///
/// If all channels use the binned likelihood optimization, the compilation is
/// delegated to compileSimPdfAsBinnedMixture() above.
///
/// Returns nullptr if some feature of this RooSimultaneous or of the fit
/// configuration is not supported yet, in which case the caller falls back to
/// the channel-splitting compilation.
std::unique_ptr<RooAbsArg>
compileSimPdfAsMixture(RooSimultaneous const &simPdf, RooArgSet const &normSet, RooFit::Detail::CompileContext &ctx)
{
   auto fallBack = [&](std::string const &why) {
      oocoutI(&simPdf, Fitting) << "RooSimultaneous::compileForNormSet(" << simPdf.GetName()
                                << "): not compiling to a mixture pdf (falling back to channel splitting): " << why
                                << std::endl;
      return std::unique_ptr<RooAbsArg>{};
   };

   RooAbsCategoryLValue const &indexCat = simPdf.indexCat();

   if (!indexCat.isFundamental() || !normSet.find(indexCat)) {
      return fallBack("the index category is not a fundamental observable");
   }
   if (simPdf.getStringAttribute("RangeName")) {
      return fallBack("ranged fits are not supported yet");
   }
   for (RooAbsArg *arg : normSet) {
      if (arg->getAttribute("__conditional__")) {
         return fallBack("conditional observables (projected dependents) are not supported yet");
      }
   }

   // Classify the channels upfront: if all of them use the binned likelihood
   // optimization, the simultaneous pdf is compiled into a single
   // concatenated binned likelihood instead of a normalized mixture pdf.
   // Mixed binned/unbinned configurations are not supported.
   std::size_t nChannels = 0;
   std::size_t nBinnedL = 0;
   for (auto const &catState : indexCat) {
      RooAbsPdf *channelPdf = simPdf.getPdf(catState.first.c_str());
      if (!channelPdf) {
         // The channel-splitting path silently drops data entries of states
         // without an associated pdf. The mixture pdf would evaluate to zero
         // for them, so it can't reproduce that behavior.
         return fallBack("state \"" + catState.first + "\" has no pdf attached");
      }
      ++nChannels;
      if (RooHelpers::getBinnedL(*channelPdf).isBinnedL) {
         ++nBinnedL;
      }
   }
   if (nChannels == 0) {
      return fallBack("there are no channels");
   }
   if (nBinnedL == nChannels) {
      return compileSimPdfAsBinnedMixture(simPdf, normSet, ctx, fallBack);
   }
   if (nBinnedL > 0) {
      return fallBack("mixed binned-likelihood and unbinned channels are not supported yet");
   }
   if (ctx.binOffsetMode()) {
      // The bin-by-bin offsetting of RooNLLVarNew builds a template pdf from
      // the dataset that is normalized over all events, while the
      // channel-splitting path uses per-channel templates with per-channel
      // weight-sum normalizations. Reproducing that in the mixture would need
      // an offset template that is conditional on the index variable. Note
      // that the all-binned compilation above doesn't have this problem: the
      // binned likelihood offsets each Poisson term with the observed bin
      // content directly, with no template pdf involved.
      return fallBack("bin-by-bin likelihood offsetting is not supported yet");
   }

   // The stand-in for the index category. The range is set once all state
   // indices are known.
   auto standIn = std::make_unique<RooRealVar>(indexCat.GetName(), indexCat.GetTitle(), 0.0);

   RooArgList prods;
   bool allExtendable = true;
   int minIndex = std::numeric_limits<int>::max();
   int maxIndex = std::numeric_limits<int>::min();

   std::map<std::string, RooAbsPdf const *> seenChannelPdfs;

   for (auto const &catState : indexCat) {
      RooAbsPdf *channelPdf = simPdf.getPdf(catState.first.c_str());
      auto [seenIt, inserted] = seenChannelPdfs.emplace(channelPdf->GetName(), channelPdf);
      if (!inserted && seenIt->second != channelPdf) {
         // Two different pdf objects with the same name would collide in
         // the name-keyed deduplication of the graph compilation. Attaching
         // the same pdf object to several channels is fine: also the
         // identically-named expected-events functions of an extended fit
         // are deduplicated to one shared clone.
         return fallBack("two channels use different pdfs with the same name \"" + seenIt->first + "\"");
      }
      allExtendable &= channelPdf->canBeExtended();
      minIndex = std::min(minIndex, catState.second);
      maxIndex = std::max(maxIndex, catState.second);

      // The suffixes make the names collision-safe: a plain "<sim>_<cat>"
      // can easily coincide with the name of an existing pdf (e.g. a channel
      // pdf named after the simultaneous pdf and the channel), and duplicate
      // names in the compiled computation graph are not allowed.
      std::string baseName = std::string(simPdf.GetName()) + "_" + catState.first;
      auto indicator = std::make_unique<RooFit::Detail::RooChannelIndicatorPdf>(
         (baseName + "_mixtureIndicator").c_str(), (baseName + "_mixtureIndicator").c_str(), *standIn, catState.second);
      // Declare to the RooFit::Evaluator that this node is a data-only
      // {0,1}-valued mask: the evaluator can then restrict the evaluation of
      // the other factors in the gated product to the events selected by the
      // mask (see Evaluator::rangeRestrictionAnalysis()). The attribute is
      // copied along when the node is cloned during graph compilation.
      indicator->setAttribute("BinaryMask");
      auto prod = std::make_unique<RooProdPdf>((baseName + "_mixtureTerm").c_str(), (baseName + "_mixtureTerm").c_str(),
                                               RooArgList(*indicator, *channelPdf));
      prod->addOwnedComponents(std::move(indicator));
      prods.addOwned(std::move(prod));
   }

   if (ctx.extendedMode() && !allExtendable) {
      return fallBack("extended fit with non-extendable channel pdfs");
   }

   standIn->setRange(minIndex - 0.5, maxIndex + 0.5);
   standIn->setVal(indexCat.getCurrentIndex());

   std::unique_ptr<RooAddPdf> mixture;
   std::unique_ptr<RooConstVar> coefVar;
   if (ctx.extendedMode()) {
      // All-extendable mode: the coefficients are the expected event yields.
      mixture = std::make_unique<RooAddPdf>(simPdf.GetName(), simPdf.GetTitle(), prods);
   } else {
      // Constant coefficients 1/C, matching the legacy convention of adding
      // sumOfWeights * log(nChannels) to the simultaneous NLL exactly. An
      // owned constant is used instead of RooFit::RooConst(), because the
      // global constants registry must not end up in a compiled computation
      // graph (concurrent evaluators would clash on its data token).
      std::string coefName = std::string(simPdf.GetName()) + "_mixtureCoef";
      coefVar = std::make_unique<RooConstVar>(coefName.c_str(), coefName.c_str(), 1.0 / prods.size());
      RooArgList coefs;
      for (std::size_t i = 0; i + 1 < prods.size(); ++i) {
         coefs.add(*coefVar);
      }
      mixture = std::make_unique<RooAddPdf>(simPdf.GetName(), simPdf.GetTitle(), prods, coefs);
   }

   // The normalization set for the mixture, with the index category replaced
   // by its real-valued stand-in.
   RooArgSet mixtureNormSet;
   for (RooAbsArg *arg : normSet) {
      mixtureNormSet.add(arg->namePtr() == indexCat.namePtr() ? *standIn : *arg);
   }

   mixture->addOwnedComponents(std::move(standIn));
   mixture->addOwnedComponents(std::move(prods));
   if (coefVar) {
      mixture->addOwnedComponents(std::move(coefVar));
   }

   std::unique_ptr<RooAbsArg> compiled = mixture->compileForNormSet(mixtureNormSet, ctx);

   if (ctx.extendedMode()) {
      // In the non-extended case, the constant mixture weights 1/C reproduce
      // the legacy convention of adding sumOfWeights * log(nChannels) to the
      // simultaneous NLL. In the extended case, the mixture weights are the
      // expected yields instead, so the term has to be requested from the
      // likelihood class explicitly to get exactly the same NLL values.
      compiled->setStringAttribute("SimCount", std::to_string(nChannels).c_str());
   }

   // Mark the compiled mixture terms as products gated by a binary mask, so
   // that the RooFit::Evaluator can restrict the evaluation of the channel
   // pdfs to the events of their own channel. The compiled product nodes are
   // new objects that don't inherit the attributes of the RooProdPdfs above,
   // so the marking has to happen after the compilation.
   if (auto *compiledAddPdf = dynamic_cast<RooAddPdf *>(compiled.get())) {
      for (RooAbsArg *component : compiledAddPdf->pdfList()) {
         for (RooAbsArg *server : component->servers()) {
            if (server->getAttribute("BinaryMask") && server->isValueServer(*component)) {
               component->setAttribute("MaskGatedProduct");
               break;
            }
         }
      }
   }

   // Keep the uncompiled mixture template alive: some normalization sets
   // stored inside RooProdPdf are disconnected from the computation graph, so
   // server redirection has no control over them (see the comment in the
   // channel-splitting compilation below). Rename it to avoid a name clash
   // with the compiled pdf.
   mixture->SetName((std::string("_") + mixture->GetName()).c_str());
   compiled->addOwnedComponents(std::move(mixture));

   return compiled;
}

void markObs(RooAbsArg *arg, std::string const &prefix, RooArgSet const &normSet)
{
   for (RooAbsArg *server : arg->servers()) {
      if (server->isFundamental() && normSet.find(*server)) {
         markObs(server, prefix, normSet);
         server->setAttribute("__obs__");
      } else if (!server->isFundamental()) {
         markObs(server, prefix, normSet);
      }
   }
}

void prefixArgs(RooAbsArg *arg, std::string const &prefix, RooArgSet const &normSet)
{
   if (!arg->getStringAttribute("__prefix__")) {
      arg->SetName((prefix + arg->GetName()).c_str());
      arg->setStringAttribute("__prefix__", prefix.c_str());
   }
   for (RooAbsArg *server : arg->servers()) {
      if (server->isFundamental() && normSet.find(*server)) {
         prefixArgs(server, prefix, normSet);
      } else if (!server->isFundamental()) {
         prefixArgs(server, prefix, normSet);
      }
   }
}

} // namespace

std::unique_ptr<RooAbsArg>
RooSimultaneous::compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext &ctx) const
{
   if (!indexCatIsObservable(normSet)) {
      // The index category is not an observable here: the RooSimultaneous
      // acts as a plain "switch" that evaluates to the component selected by
      // the current index state, analogous to RooMultiPdf. The channel
      // observables are then the same as the ones of this pdf, so the
      // channel-splitting compilation below (which renames the per-channel
      // observables so they can be filled from split datasets) must not be
      // used. Compile like an ordinary self-normalized pdf instead.
      return RooAbsPdf::compileForNormSet(normSet, ctx);
   }

   // Experimental, opt-in via ROOFIT_SIM_COMPILE_MIXTURE=1: compile into an
   // ordinary mixture pdf built from standard components, so that the
   // downstream likelihood machinery needs no special treatment of the
   // simultaneous case. Only for likelihood compilation; unsupported
   // configurations fall through to the channel-splitting path below.
   if (ctx.likelihoodMode() && simMixtureCompileRequested()) {
      if (std::unique_ptr<RooAbsArg> mixture = compileSimPdfAsMixture(*this, normSet, ctx)) {
         return mixture;
      }
   }

   std::unique_ptr<RooSimultaneous> newSimPdf{static_cast<RooSimultaneous *>(this->Clone())};

   const char *rangeName = this->getStringAttribute("RangeName");
   bool splitRange = this->getAttribute("SplitRange");

   RooArgSet newPdfs;
   std::vector<std::string> catNames;

   for (auto *proxy : static_range_cast<RooRealProxy *>(newSimPdf->_pdfProxyList)) {
      catNames.emplace_back(proxy->GetName());
      std::string const &catName = catNames.back();
      const std::string prefix = "_" + catName + "_";

      const std::string origname = proxy->arg().GetName();

      auto pdfClone = RooHelpers::cloneTreeWithSameParameters(static_cast<RooAbsPdf const &>(proxy->arg()), &normSet);

      markObs(pdfClone.get(), prefix, normSet);

      std::unique_ptr<RooArgSet> pdfNormSet{
         std::unique_ptr<RooArgSet>(pdfClone->getVariables())->selectByAttrib("__obs__", true)};
      std::unique_ptr<RooArgSet> condVarSet{
         std::unique_ptr<RooArgSet>(pdfClone->getVariables())->selectByAttrib("__conditional__", true)};

      pdfNormSet->remove(*condVarSet, true, true);

      if (rangeName) {
         pdfClone->setNormRange(RooHelpers::getRangeNameForSimComponent(rangeName, splitRange, catName).c_str());
      }

      RooFit::Detail::CompileContext pdfContext{*pdfNormSet};
      pdfContext.setLikelihoodMode(ctx.likelihoodMode());
      auto *pdfFinal = pdfContext.compile(*pdfClone, *newSimPdf, *pdfNormSet);

      // We can only prefix the observables after everything related the
      // compiling of the compute graph for the normalization set is done. This
      // is because of a subtlety in conditional RooProdPdfs, which stores the
      // normalization sets for the individual pdfs in RooArgSets that are
      // disconnected from the computation graph, so we have no control over
      // them. An alternative would be to use recursive server re-direction,
      // but this has more performance overhead.
      prefixArgs(pdfFinal, prefix, normSet);

      pdfFinal->fixAddCoefNormalization(*pdfNormSet, false);

      pdfClone->SetName((std::string("_") + pdfClone->GetName()).c_str());
      pdfFinal->addOwnedComponents(std::move(pdfClone));

      pdfFinal->setAttribute(("ORIGNAME:" + origname).c_str());
      newPdfs.add(*pdfFinal);

      // We will remove the old pdf server because we will fill the new ones by
      // hand via the creation of new proxies.
      newSimPdf->removeServer(const_cast<RooAbsReal &>(proxy->arg()), true);
   }

   // Replace pdfs with compiled pdfs. Don't use RooAbsArg::redirectServers()
   // here, because it doesn't support replacing two servers with the same name
   // (it can happen in a RooSimultaneous that two pdfs have the same name).

   // First delete old proxies (we have already removed the servers before).
   newSimPdf->_pdfProxyList.Delete();

   // Recreate the _pdfProxyList with the compiled pdfs
   for (std::size_t i = 0; i < newPdfs.size(); ++i) {
      const char *label = catNames[i].c_str();
      newSimPdf->_pdfProxyList.Add(
         new RooRealProxy(label, label, newSimPdf.get(), *static_cast<RooAbsReal *>(newPdfs[i])));
   }

   ctx.compileServers(*newSimPdf, normSet); // to trigger compiling also the index category

   return newSimPdf;
}
