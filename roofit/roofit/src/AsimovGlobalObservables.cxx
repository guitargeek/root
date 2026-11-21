/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2026
 *
 * Copyright (c) 2026, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <RooFit/AsimovGlobalObservables.h>

#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooBifurGauss.h>
#include <RooGamma.h>
#include <RooGaussian.h>
#include <RooLognormal.h>
#include <RooMsgService.h>
#include <RooPoisson.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>
#include <RooUniform.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

const char *const prefix = "RooFit::setGlobalObservablesToAsimov: constraint term ";

// Of the two arguments of a constraint term that it relates to each other, like the x and the mean of a Gaussian, the
// one that is not the global observable. A `nullptr` if none of them is the global observable.
RooAbsReal const *otherSide(RooAbsReal const &a, RooAbsReal const &b, RooRealVar const &glob)
{
   if (&a == &glob) {
      return &b;
   }
   if (&b == &glob) {
      return &a;
   }
   return nullptr;
}

// For constraint terms that are not known by their class: the server of the term that stands for the constrained
// parameter. It is the only server that doesn't depend on a global observable and has a floating parameter. If the
// parameter is constant, it is the only such server that has parameters at all.
RooAbsReal const *parameterSide(RooAbsPdf const &constraint, RooArgSet const &globalObservables, bool &ambiguous)
{
   std::vector<RooAbsReal const *> floating;
   std::vector<RooAbsReal const *> fixed;
   for (RooAbsArg const *server : constraint.servers()) {
      auto const *real = dynamic_cast<RooAbsReal const *>(server);
      if (!real || real->dependsOn(globalObservables)) {
         continue;
      }
      // The parameters of the server: itself, if it is a variable.
      RooArgSet parameters;
      if (real->isFundamental()) {
         parameters.add(*real);
      } else {
         real->getParameters(nullptr, parameters);
      }
      if (parameters.empty()) {
         continue; // a number
      }
      const bool isFloating =
         std::any_of(parameters.begin(), parameters.end(), [](RooAbsArg const *arg) { return !arg->isConstant(); });
      (isFloating ? floating : fixed).push_back(real);
   }
   ambiguous = floating.size() > 1 || (floating.empty() && fixed.size() > 1);
   if (floating.size() == 1) {
      return floating[0];
   }
   return floating.empty() && fixed.size() == 1 ? fixed[0] : nullptr;
}

void setToAsimov(RooAbsPdf &constraint, RooArgSet const &globalObservables)
{
   // A product that groups constraint terms.
   if (auto *prod = dynamic_cast<RooProdPdf *>(&constraint)) {
      for (auto *factor : static_range_cast<RooAbsPdf *>(prod->pdfList())) {
         setToAsimov(*factor, globalObservables);
      }
      return;
   }
   // A uniform term constrains nothing.
   if (dynamic_cast<RooUniform *>(&constraint)) {
      return;
   }

   RooArgSet globs;
   constraint.getObservables(&globalObservables, globs);
   if (globs.empty()) {
      oocoutW(nullptr, Generation) << prefix << constraint.GetName() << " has no global observables - skip it"
                                   << std::endl;
      return;
   }
   auto *glob = dynamic_cast<RooRealVar *>(globs.first());
   if (globs.size() > 1 || !glob) {
      oocoutE(nullptr, Generation) << prefix << constraint.GetName()
                                   << " has multiple global observables - cannot set them - skip it" << std::endl;
      return;
   }

   // The terms that relate two of their arguments to each other: the global observable takes the value of the other
   // one.
   RooAbsReal const *other = nullptr;
   bool isKnown = true;
   bool isLogOfOther = false;
   bool isExpOfOther = false;
   if (auto *gauss = dynamic_cast<RooGaussian *>(&constraint)) {
      other = otherSide(gauss->getX(), gauss->getMean(), *glob);
   } else if (auto *bifur = dynamic_cast<RooBifurGauss *>(&constraint)) {
      other = otherSide(bifur->getX(), bifur->getMean(), *glob);
   } else if (auto *pois = dynamic_cast<RooPoisson *>(&constraint)) {
      // The global observable takes the value of the mean, which is no integer in general.
      pois->setNoRounding(true);
      other = otherSide(pois->getX(), pois->getMean(), *glob);
   } else if (auto *logn = dynamic_cast<RooLognormal *>(&constraint)) {
      other = otherSide(logn->getX(), logn->getMedian(), *glob);
      // In the standard parametrization, the parameter is the logarithm of the median.
      isExpOfOther = logn->useStandardParametrization() && other == &logn->getMedian();
      isLogOfOther = logn->useStandardParametrization() && other == &logn->getX();
   } else {
      isKnown = false;
   }
   if (other) {
      const double value = other->getVal();
      glob->setVal(isExpOfOther ? std::exp(value) : (isLogOfOther ? std::log(value) : value));
      return;
   }
   if (isKnown) {
      oocoutE(nullptr, Generation) << prefix << constraint.GetName()
                                   << " has no direct dependence on global observable - cannot set it" << std::endl;
      return;
   }

   // Other terms: the global observable takes the value of the server that stands for the parameter. This is
   // right for the Gamma constraint of HistFactory, Gamma(b | k = n + 1, theta), whose global observable n is the
   // mode of b / theta. The global observable is no direct server there, as k is a function of it.
   auto *gamma = dynamic_cast<RooGamma *>(&constraint);
   if (!gamma) {
      oocoutW(nullptr, Generation) << prefix << constraint.GetName() << " of type " << constraint.ClassName()
                                   << " is a non-supported type - result might be not correct" << std::endl;
      if (!constraint.findServer(*glob)) {
         oocoutE(nullptr, Generation) << prefix << constraint.GetName()
                                      << " has no direct dependence on global observable - cannot set it" << std::endl;
         return;
      }
   }
   bool ambiguous = false;
   RooAbsReal const *parameter = parameterSide(constraint, globalObservables, ambiguous);
   if (!parameter) {
      oocoutE(nullptr, Generation) << prefix << constraint.GetName()
                                   << (ambiguous ? " has more than one server that depends on parameters"
                                                 : " has no server that depends on parameters")
                                   << " - the global observable " << glob->GetName()
                                   << " will not be set to its Asimov value" << std::endl;
      return;
   }
   const double scale = gamma ? gamma->getBeta().getVal() : 1.0;
   glob->setVal(scale > 0. ? parameter->getVal() / scale : parameter->getVal());
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Set the global observables of constraint terms to their Asimov values: the
/// values that the auxiliary measurements are expected to give for the current
/// values of the parameters, such that the constraint terms pull none of the
/// parameters away from them.
///
/// For a RooGaussian, RooBifurGauss, RooPoisson or RooLognormal that relates a
/// global observable to a parameter (or a function of parameters), the global
/// observable takes the value of the parameter: these terms have the global
/// observable as `x` and the parameter as mean or median, or the other way
/// round. For a RooLognormal in the standard parametrization, where the
/// parameter is the logarithm of the median, it takes the corresponding value.
/// The rounding of a RooPoisson is switched off, see
/// RooPoisson::setNoRounding(), as the value is no integer in general.
///
/// For a RooGamma, the global observable takes the value of the parameter
/// divided by the scale `beta`, which is right for the Gamma constraints of
/// HistFactory, Gamma(b | k = n + 1, theta) with the global observable n. For
/// terms of other types, it takes the value of the only argument of the term
/// that depends on parameters, with a warning, as this may not be what the
/// term means.
///
/// Terms that are products of constraint terms are set through their factors.
/// Terms with no or with several global observables are skipped with a
/// message.
///
/// \param[in] constraints The constraint terms, e.g. the external constraints
///            of a ModelConfig.
/// \param[in] globalObservables The global observables of the model. The ones
///            that the constraint terms depend on are changed, also if they are
///            constant.
///
/// \see RooFit::generateGlobalObservables() to sample the global observables
///      instead.

void RooFit::setGlobalObservablesToAsimov(RooAbsCollection const &constraints, RooArgSet const &globalObservables)
{
   for (RooAbsArg *term : constraints) {
      if (auto *constraint = dynamic_cast<RooAbsPdf *>(term)) {
         setToAsimov(*constraint, globalObservables);
      }
   }
}
