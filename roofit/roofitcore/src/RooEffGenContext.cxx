/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven, NIKHEF/VU, Gerhard.Raven@nikhf.nl                    *
 *                                                                           *
 * Copyright (c) 2005, NIKHEF.  All rights reserved.                         *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooEffGenContext.cxx
\class RooEffGenContext
\ingroup Roofitcore

Specialized generator context for p.d.fs represented
by class RooEffProd, which are p.d.fs multiplied with an efficiency function.
This generator context generates events from such products by first
generating events from a dedicated generator context of the input p.d.f.
and applying an extra rejection step based on the efficiency function.
**/

#include <memory>

#include "RooEffGenContext.h"
#include "RooAbsPdf.h"
#include "RooRandom.h"


////////////////////////////////////////////////////////////////////////////////
/// Constructor of generator context for RooEffProd products

RooEffGenContext::RooEffGenContext(const RooAbsPdf &model, const RooAbsPdf &pdf, const RooAbsReal &eff,
                                   const RooArgSet &vars, const RooDataSet *prototype, const RooArgSet *auxProto,
                                   bool verbose, const RooArgSet * /*forceDirect*/)
   : RooAbsGenContext(model, vars, prototype, auxProto, verbose), _maxEff(0.)
{
   RooArgSet x(eff, eff.GetName());
   x.snapshot(_cloneSet, true);
   initializeEff(eff);
   _generator = std::unique_ptr<RooAbsGenContext>{pdf.genContext(vars, prototype, auxProto, verbose)};
   vars.snapshot(_vars, true);
}

////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of generator.

void RooEffGenContext::initGenerator(const RooArgSet &theEvent)
{
   _eff->recursiveRedirectServers(theEvent);
   _generator->initGenerator(theEvent);

   // Check if PDF supports maximum finding
   Int_t code = _eff->getMaxVal(_vars);
   if (!code) {
      _maxEff = 1.;
   } else {
      _maxEff = _eff->maxVal(code);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generate one event. Generate an event from the p.d.f and
/// then perform an accept/reject sampling based on the efficiency
/// function
///
/// \note This accept/reject step still handles a violated maximum badly, and
/// this was left alone on purpose to keep the fix for issue #12317 focused.
/// If the efficiency does not advertise a maximum via
/// RooAbsReal::getMaxVal(), a maximum of one is assumed, and any point where
/// the efficiency exceeds that is rejected outright instead of the maximum
/// being adjusted, which zeroes out exactly the region with the highest
/// efficiency (and does not even terminate if the efficiency does not depend
/// on the generated observables). If a maximum *is* advertised but is too
/// small, every point above it is accepted, which truncates the distribution
/// like it used to in RooAcceptReject. Both cases should be fixed the way
/// RooAcceptReject does it, by raising the maximum and resampling what was
/// already generated.

void RooEffGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
   while (true) {
      _generator->generateEvent(theEvent, remaining);
      double val = _eff->getVal();
      if (val > _maxEff && !_eff->getMaxVal(_vars)) {
         coutE(Generation) << ClassName() << "::" << GetName()
                           << ":generateEvent: value of efficiency is larger than assumed maximum of 1." << std::endl;
         continue;
      }
      if (val > RooRandom::uniform() * _maxEff) {
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooEffGenContext::printMultiline(std::ostream &os, Int_t content, bool verbose, TString indent) const
{
   RooAbsGenContext::printMultiline(os, content, verbose, indent);
   os << indent << "--- RooEffGenContext ---" << std::endl;
   os << indent << "Using EFF ";
   _eff->printStream(os, kName | kArgs | kClassName, kSingleLine, indent);
   os << indent << "PDF generator" << std::endl;

   TString indent2(indent);
   indent2.Append("    ");

   _generator->printMultiline(os, content, verbose, indent2);
}
