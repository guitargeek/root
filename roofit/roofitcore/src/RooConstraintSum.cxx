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
\file RooConstraintSum.cxx
\class RooConstraintSum
\ingroup Roofitcore

Calculates the sum of the -(log) likelihoods of
a set of RooAbsPfs that represent constraint functions. This class
is used to calculate the composite -log(L) of constraints to be
added to the regular -log(L) in RooAbsPdf::fitTo() with Constrain(..)
arguments.
**/


#include "RooConstraintSum.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooHelpers.h"
#include "RooAbsCategoryLValue.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooRealVar.h"



////////////////////////////////////////////////////////////////////////////////
/// Constructor with set of constraint p.d.f.s. All elements in constraintSet must inherit from RooAbsPdf.

RooConstraintSum::RooConstraintSum(const char* name, const char* title, const RooArgSet& constraintSet, const RooArgSet& normSet, bool takeGlobalObservablesFromData) :
  RooAbsReal(name, title),
  _gaussianParams("gaussianParams","gaussianParams",this),
  _gaussianObs("gaussianObs","gaussianObs",this),
  _poissonParams("poissonParams","poissonParams",this),
  _poissonObs("poissonObs","poissonObs",this),
  _set1("set1","First set of components",this),
  _takeGlobalObservablesFromData{takeGlobalObservablesFromData}
{
  _set1.addTyped<RooAbsPdf>(constraintSet);
  _paramSet.add(normSet) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

RooConstraintSum::RooConstraintSum(const RooConstraintSum& other, const char* name) :
  RooAbsReal(other, name),
  _gaussians{other._gaussians},
  _gaussianParams("gaussianParams",this,other._gaussianParams),
  _gaussianObs("gaussianObs",this,other._gaussianObs),
  _poissons{other._poissons},
  _poissonParams("poissonParams",this,other._poissonParams),
  _poissonObs("poissonObs",this,other._poissonObs),
  _set1("set1",this,other._set1),
  _paramSet(other._paramSet),
  _takeGlobalObservablesFromData{other._takeGlobalObservablesFromData}
{
}


////////////////////////////////////////////////////////////////////////////////
/// Return sum of -log of constraint p.d.f.s.

double RooConstraintSum::evaluate() const
{
  double sum(0);

  for (const auto comp : _set1) {
    sum -= static_cast<RooAbsPdf*>(comp)->getLogVal(&_paramSet);
  }

  for (std::size_t iGauss = 0; iGauss < _gaussians.size(); ++iGauss) {
    auto const &info = _gaussians[iGauss];
    const double paramVal = static_cast<RooAbsReal const &>(_gaussianParams[iGauss]).getVal();
    const double obsVal = static_cast<RooAbsReal const &>(_gaussianObs[iGauss]).getVal();
    const double arg = info.sigmaInvVal * (obsVal - paramVal);
    const double lnVal = std::log(info.sigmaInvVal / std::sqrt(2 * M_PI)) - 0.5 * arg * arg;
    sum -= lnVal;
  }

  for (std::size_t iPoiss = 0; iPoiss < _poissons.size(); ++iPoiss) {
    auto const &info = _poissons[iPoiss];
    const double paramVal = static_cast<RooAbsReal const &>(_poissonParams[iPoiss]).getVal();
    const double obsVal = static_cast<RooAbsReal const &>(_poissonObs[iPoiss]).getVal();
    double k = info.noRounding ? obsVal : std::floor(obsVal);
    sum -= k * std::log(paramVal) - std::lgamma(k + 1.) - paramVal;
  }

  return sum;
}

void RooConstraintSum::doEval(RooFit::EvalContext &ctx) const
{
   double sum(0);

   for (const auto comp : _set1) {
      sum -= std::log(ctx.at(comp)[0]);
   }

   for (std::size_t iGauss = 0; iGauss < _gaussians.size(); ++iGauss) {
      auto const &info = _gaussians[iGauss];
      const double paramVal = ctx.at(&_gaussianParams[iGauss])[0];
      const double obsVal = ctx.at(&_gaussianObs[iGauss])[0];
      const double arg = info.sigmaInvVal * (obsVal - paramVal);
      const double lnVal = std::log(info.sigmaInvVal / std::sqrt(2 * M_PI)) - 0.5 * arg * arg;
      sum -= lnVal;
   }

   for (std::size_t iPoiss = 0; iPoiss < _poissons.size(); ++iPoiss) {
      auto const &info = _poissons[iPoiss];
      const double paramVal = ctx.at(&_poissonParams[iPoiss])[0];
      const double obsVal = ctx.at(&_poissonObs[iPoiss])[0];
      double k = info.noRounding ? obsVal : std::floor(obsVal);
      sum -= k * std::log(paramVal) - std::lgamma(k + 1.) - paramVal;
   }

   ctx.output()[0] = sum;
}

bool RooConstraintSum::hardcodeGaussian(RooGaussian const &gauss)
{
   const bool xIsObs = _paramSet.contains(gauss.getX());
   const bool meanIsObs = _paramSet.contains(gauss.getMean());

   // The actual mean parameter and observable in this context
   RooRealVar const *obs = dynamic_cast<RooRealVar const *>(xIsObs ? &gauss.getX() : &gauss.getMean());
   RooRealVar const *mu = dynamic_cast<RooRealVar const *>(xIsObs ? &gauss.getMean() : &gauss.getX());

   // Don't attempt to hardcode the cases where both "x" and "mean" are
   // observables, or where the parameter or observable is not a RooRealVar.
   if (xIsObs == meanIsObs || !mu || !obs)
      return false;

   // The width value is baked into the cached info, so it must not change
   // during the lifetime of this object.
   if (!gauss.getSigma().isConstant())
      return false;

   GaussianInfo info;
   info.sigmaInvVal = 1. / gauss.getSigma().getVal();

   const double howManySigmas = std::min(mu->getVal() - mu->getMin(), mu->getMax() - mu->getVal()) * info.sigmaInvVal;

   // A good threshold would be 5*sigma
   if (howManySigmas < 4.99)
      return false;

   _gaussians.emplace_back(info);
   _gaussianParams.add(*mu);
   _gaussianObs.add(*obs);
   return true;
}

bool RooConstraintSum::hardcodePoisson(RooPoisson const &poiss)
{
   PoissonInfo info;
   info.noRounding = poiss.getNoRounding();
   _poissons.emplace_back(info);
   _poissonParams.add(poiss.getMean());
   _poissonObs.add(poiss.getX());
   return true;
}

std::unique_ptr<RooAbsArg> RooConstraintSum::compileForNormSet(RooArgSet const & /*normSet*/, RooFit::Detail::CompileContext & ctx) const
{
   auto compiledConstraints =
      std::make_unique<RooConstraintSum>(GetName(), GetTitle(), RooArgSet{}, _paramSet, _takeGlobalObservablesFromData);

   for (RooAbsArg *arg : _set1) {
      bool hardcoded = false;
      if (auto gauss = dynamic_cast<RooGaussian *>(arg)) {
         hardcoded = compiledConstraints->hardcodeGaussian(*gauss);
      } else if (auto poiss = dynamic_cast<RooPoisson *>(arg)) {
         hardcoded = compiledConstraints->hardcodePoisson(*poiss);
      }

      if (!hardcoded) {
         compiledConstraints->_set1.add(*arg);
      }
   }

   for (const auto server : compiledConstraints->servers()) {
      RooArgSet nset;
      server->getObservables(&_paramSet, nset);
      ctx.compileServer(*server, *compiledConstraints, nset);
   }

   return compiledConstraints;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace the variables in this RooConstraintSum with the global observables
/// in the dataset if they match by name. This function will do nothing if this
/// RooConstraintSum is configured to not use the global observables stored in
/// datasets.
bool RooConstraintSum::setData(RooAbsData const& data, bool /*cloneData=true*/) {
  if(_takeGlobalObservablesFromData && data.getGlobalObservables()) {
    this->recursiveRedirectServers(*data.getGlobalObservables()) ;
  }
  return true;
}
