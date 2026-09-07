/// \cond ROOFIT_INTERNAL

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl   *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
/// \class RooMinimizerFcn
/// RooMinimizerFcn is an interface to the ROOT::Math::IBaseFunctionMultiDim,
/// a function that ROOT's minimisers use to carry out minimisations.
///

#include "RooMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooArgSet.h"
#include "RooConstraintSum.h"
#include "RooEvaluatorWrapper.h"
#include "RooFit/Detail/RooChannelIndicatorPdf.h"
#include "RooFit/Detail/RooNLLVarNew.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooNaNPacker.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMatrixDSym.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <vector>

using std::setprecision;

namespace {

// Check whether two sorted ranges have at least one element in common.
// Like std::set_intersection, both input ranges must be sorted; the early
// return on the first match keeps the common case cheap.
template <class InputIt1, class InputIt2>
bool intersect(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
{
   while (first1 != last1 && first2 != last2) {
      if (*first1 < *first2) {
         ++first1;
         continue;
      }
      if (*first2 < *first1) {
         ++first2;
         continue;
      }
      return true;
   }
   return false;
}

// Helper function that wraps RooAbsArg::getParameters and directly returns the
// output RooArgSet. To be used in the initializer list of the RooMinimizerFcn
// constructor. In the case of figuring out all parameters for the minimizer,
// we don't want to strip disconnected parameters, becuase which parameters are
// disconnected can change between minimization runs.
RooArgSet getAllParameters(RooAbsReal const &funct)
{
   RooArgSet out;
   funct.getParameters(nullptr, out, /*stripDisconnected*/ false);
   return out;
}

// Groups of computation-graph leaves by the additive term of the minimized
// function they appear in. Two parameters that share no term index have a
// mixed second derivative that is identically zero.
//
// Note that distinct objects with the same name share their pointer from the
// name registry, so they end up merged in the same map entry. This is
// conservative: it can only add co-occurrences, never remove any.
struct VariableGroups {
   /// For each leaf (keyed by its unique name pointer), the sorted list of
   /// indices of the additive terms it appears in.
   std::unordered_map<TNamed const *, std::vector<int>> groups;

   /// Register one additive term: record for every leaf in the collection
   /// that it appears in this term.
   void registerTerm(RooAbsCollection const &leaves)
   {
      for (RooAbsArg const *arg : leaves) {
         groups[arg->namePtr()].push_back(_nextIndex);
      }
      ++_nextIndex;
   }

private:
   int _nextIndex = 0;
};

// Add the leaves of the computation graph under `arg` to `out`. The detour
// via RooArgList avoids the deduplication done after adding each element.
void addLeaves(RooAbsArg const &arg, RooArgSet &out)
{
   RooArgList leafList;
   arg.treeNodeServerList(&leafList, nullptr, /*branches*/ false, /*leaves*/ true, /*valueOnly*/ false,
                          /*recurseFundamental*/ true);
   out.add(leafList.begin(), leafList.end());
}

// If `arg` is the mixture that RooSimultaneous::compileForNormSet() compiles a
// simultaneous pdf into, return it, and nullptr otherwise.
//
// The compiled form is an ordinary mixture over the channels,
// \f$ \sum_s c_s \, \mathbf{1}[c = s] \, f_s(x) \f$, i.e. a RooRealSumPdf
// whose every func is gated by a RooChannelIndicatorPdf.
RooRealSumPdf const *channelMixturePdf(RooAbsArg const &arg)
{
   auto const *sumPdf = dynamic_cast<RooRealSumPdf const *>(&arg);
   if (!sumPdf || sumPdf->funcList().size() < 2) {
      return nullptr;
   }
   for (RooAbsArg const *func : sumPdf->funcList()) {
      auto isIndicator = [](RooAbsArg const *server) {
         return dynamic_cast<RooFit::Detail::RooChannelIndicatorPdf const *>(server) != nullptr;
      };
      if (std::none_of(func->servers().begin(), func->servers().end(), isIndicator)) {
         return nullptr;
      }
   }
   return sumPdf;
}

// Register one additive term per channel for a likelihood over a
// mixture-compiled simultaneous pdf, and report whether `arg` was such a
// likelihood.
//
// The channel indicators are mutually exclusive: for any index state exactly
// one of them is one and all others are zero. Every data row therefore
// contributes to exactly one channel, and the likelihood is the strictly
// additive sum over the channels of their per-channel likelihoods. This is
// the mixture-compiled equivalent of the sum over per-channel terms that a
// RooAddition exposes directly.
bool fillMixtureVariableGroups(RooAbsArg const &arg, VariableGroups &out)
{
   if (!dynamic_cast<RooFit::Detail::RooNLLVarNew const *>(&arg)) {
      return false;
   }

   RooRealSumPdf const *mixture = nullptr;
   for (RooAbsArg const *server : arg.servers()) {
      if ((mixture = channelMixturePdf(*server))) {
         break;
      }
   }
   if (!mixture) {
      return false;
   }

   RooArgList const &funcs = mixture->funcList();
   RooArgList const &coefs = mixture->coefList();

   // The leaves of a channel are those of its gated term plus those of its
   // coefficient, which carries the expected number of events in extended
   // fits. A pdf shared by several channels correctly ends up in each of
   // their groups, which keeps its parameters coupled.
   std::vector<RooArgSet> channelLeaves(funcs.size());
   RooArgSet attributed;
   for (std::size_t i = 0; i < funcs.size(); ++i) {
      addLeaves(*funcs.at(i), channelLeaves[i]);
      if (i < coefs.size()) {
         addLeaves(*coefs.at(i), channelLeaves[i]);
      }
      attributed.add(channelLeaves[i], /*silent*/ true);
   }

   // Whatever else the likelihood depends on cannot be attributed to a single
   // channel, so it goes into every group to stay coupled with everything.
   // The offsetting and expected-events nodes are additive over the channels
   // as well, but they are not worth the extra bookkeeping: their parameters
   // already appear in the channels they belong to.
   RooArgSet unattributed;
   addLeaves(arg, unattributed);
   unattributed.remove(attributed, /*silent*/ true, /*matchByNameOnly*/ true);

   for (RooArgSet &leaves : channelLeaves) {
      leaves.add(unattributed, /*silent*/ true);
      out.registerTerm(leaves);
   }
   return true;
}

// Fill the map from computation-graph leaves to the additive terms of the
// minimized function they appear in.
//
// Recursing into the components of a node is only correct if the value of
// the node is a *strictly additive* combination of them: recursing into
// anything else would wrongly advertise vanishing second derivatives and
// silently corrupt Hessian results. Any other node contributes all of its
// leaves as one single term, which advertises no independence but is always
// correct.
void fillVariableGroups(RooAbsArg const &arg, VariableGroups &out)
{
   if (auto addition = dynamic_cast<RooAddition const *>(&arg)) {
      for (RooAbsArg *component : addition->list()) {
         fillVariableGroups(*component, out);
      }
      return;
   }
   if (auto constraintSum = dynamic_cast<RooConstraintSum const *>(&arg)) {
      for (RooAbsArg *component : constraintSum->list()) {
         fillVariableGroups(*component, out);
      }
      return;
   }
   if (auto wrapper = dynamic_cast<RooFit::Experimental::RooEvaluatorWrapper const *>(&arg)) {
      fillVariableGroups(wrapper->topNode(), out);
      return;
   }
   if (fillMixtureVariableGroups(arg, out)) {
      return;
   }

   RooArgSet leafSet;
   addLeaves(arg, leafSet);
   out.registerTerm(leafSet);
}

} // namespace

// use reference wrapper for the Functor, such that the functor points to this RooMinimizerFcn by reference.
RooMinimizerFcn::RooMinimizerFcn(RooAbsReal *funct, RooMinimizer *context)
   : RooAbsMinimizerFcn(getAllParameters(*funct), context), _funct(funct)
{
   unsigned int nDim = getNDim();

   if (context->_cfg.useGradient && funct->hasGradient()) {
      _gradientOutput.resize(_allParams.size());
      _multiGenFcn = std::make_unique<ROOT::Math::GradFunctor>(this, &RooMinimizerFcn::operator(),
                                                               &RooMinimizerFcn::evaluateGradient, nDim);
   } else {
      _multiGenFcn = std::make_unique<ROOT::Math::Functor>(std::cref(*this), nDim);
   }
   if (context->_cfg.useHessian) {
      _hessianOutput.resize(_allParams.size() * _allParams.size());
   }
}

/// Evaluate function given the parameters in `x`.
double RooMinimizerFcn::operator()(const double *x) const
{
   // Set the parameter values for this iteration
   for (unsigned index = 0; index < getNDim(); index++) {
      if (_logfile)
         (*_logfile) << x[index] << " ";
      SetPdfParamVal(index, x[index]);
   }

   // Calculate the function for these parameters
   RooAbsReal::setHideOffset(false);
   double fvalue = _funct->getVal();
   RooAbsReal::setHideOffset(true);

   fvalue = applyEvalErrorHandling(fvalue);

   // Optional logging
   if (_logfile)
      (*_logfile) << setprecision(15) << fvalue << setprecision(4) << std::endl;
   if (cfg().verbose) {
      std::cout << "\nprevFCN" << (_funct->isOffsetting() ? "-offset" : "") << " = " << setprecision(10) << fvalue
                << setprecision(4) << "  ";
      std::cout.flush();
   }

   finishDoEval();

   return fvalue;
}

void RooMinimizerFcn::evaluateGradient(const double *x, double *out) const
{
   // Set the parameter values for this iteration
   for (unsigned index = 0; index < getNDim(); index++) {
      if (_logfile)
         (*_logfile) << x[index] << " ";
      SetPdfParamVal(index, x[index]);
   }

   _funct->gradient(_gradientOutput.data());

   std::size_t iAll = 0;
   std::size_t iFloating = 0;
   for (RooAbsArg *param : _allParamsInit) {
      if (!treatAsConstant(*param)) {
         out[iFloating] = _gradientOutput[iAll];
         ++iFloating;
      }
      ++iAll;
   }

   // Optional logging
   if (cfg().verbose) {
      std::cout << "\n    gradient = ";
      for (std::size_t i = 0; i < getNDim(); ++i) {
         std::cout << out[i] << ", ";
      }
   }
}

std::string RooMinimizerFcn::getFunctionName() const
{
   return _funct->GetName();
}

std::string RooMinimizerFcn::getFunctionTitle() const
{
   return _funct->GetTitle();
}

void RooMinimizerFcn::setOffsetting(bool flag)
{
   _funct->enableOffsetting(flag);
}

RooArgSet RooMinimizerFcn::freezeDisconnectedParameters() const
{

   RooArgSet paramsDisconnected;
   RooArgSet paramsConnected;

   _funct->getParameters(nullptr, paramsDisconnected, /*stripDisconnected*/ false);
   _funct->getParameters(nullptr, paramsConnected, /*stripDisconnected*/ true);

   paramsDisconnected.remove(paramsConnected, true, true);

   RooArgSet changedSet;

   for (RooAbsArg *a : paramsDisconnected) {
      auto *v = dynamic_cast<RooRealVar *>(a);
      auto *cv = dynamic_cast<RooCategory *>(a);
      if (v && !v->isConstant()) {
         v->setConstant();
         changedSet.add(*v);
      } else if (cv && !cv->isConstant()) {
         cv->setConstant();
         changedSet.add(*cv);
      }
   }

   return changedSet;
}

bool RooMinimizerFcn::evaluateHessian(std::span<const double> x, double *out) const
{
   // Set the parameter values for this iteration
   for (unsigned index = 0; index < getNDim(); index++) {
      if (_logfile)
         (*_logfile) << x[index] << " ";
      SetPdfParamVal(index, x[index]);
   }

   _funct->hessian(_hessianOutput.data());

   std::size_t m = _allParamsInit.size();
   std::size_t n = getNDim();
   std::size_t iAll = 0;
   std::size_t iFloating = 0;
   for (RooAbsArg *param_i : _allParamsInit) {
      if (!treatAsConstant(*param_i)) {
         std::size_t jAll = 0;
         std::size_t jFloating = 0;
         for (RooAbsArg *param_j : _allParamsInit) {
            if (!treatAsConstant(*param_j)) {
               out[iFloating * n + jFloating] = _hessianOutput[iAll * m + jAll];
               ++jFloating;
            }
            ++jAll;
         }
         ++iFloating;
      }
      ++iAll;
   }

   // Optional logging
   if (cfg().verbose) {
      std::cout << "\n    hessian = " << std::endl;
      for (std::size_t i = 0; i < getNDim(); ++i) {
         for (std::size_t j = 0; j < getNDim(); ++j) {
            std::cout << out[i * n + j] << ", ";
         }
         std::cout << std::endl;
      }
   }
   return true;
}

void RooMinimizerFcn::initMinimizer(ROOT::Math::Minimizer &minim, RooMinimizer *context)
{
   minim.SetFunction(*_multiGenFcn);
   if (context->_cfg.useHessian && _funct->hasHessian()) {
      minim.SetHessianFunction(
         std::bind(&RooMinimizerFcn::evaluateHessian, this, std::placeholders::_1, std::placeholders::_2));
   }
   // The independence information for skipping vanishing second derivatives
   // in numerical Hessian computations is a Minuit2-only feature, so it is
   // wired up directly with the concrete minimizer type instead of going
   // through the ROOT::Math::Minimizer interface.
   if (auto *minuit2 = dynamic_cast<ROOT::Minuit2::Minuit2Minimizer *>(&minim)) {
      minuit2->SetSecondDerivativeAlwaysVanishesFunc(
         [this](unsigned int i, unsigned int j) { return secondDerivativeAlwaysVanishes(i, j); });
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the bitvector that flags for each pair of floatable parameters
/// whether they appear together in at least one additive term of the
/// minimized function, i.e. whether their mixed second derivative can be
/// non-vanishing. Built lazily because it is only needed for Hessian
/// evaluations, and building it for models with many parameters is not free.
void RooMinimizerFcn::buildSecondDerivMask() const
{
   VariableGroups groups;
   fillVariableGroups(*_funct, groups);

   std::size_t nParams = getNDim();

   // Packed bitvector: bit set means the parameter pair shares an additive
   // term, so the mixed second derivative can be non-zero.
   _secondDerivMask.assign(nParams * nParams, false);
   for (std::size_t i = 0; i < nParams; ++i) {
      _secondDerivMask[nParams * i + i] = true;
      auto found1 = groups.groups.find(floatableParam(i).namePtr());
      for (std::size_t j = 0; j < i; ++j) {
         auto found2 = groups.groups.find(floatableParam(j).namePtr());
         // A parameter that was not seen in the computation graph traversal
         // is conservatively treated as intersecting with everything.
         bool canBeNonZero = found1 == groups.groups.end() || found2 == groups.groups.end() ||
                             intersect(found1->second.begin(), found1->second.end(), found2->second.begin(),
                                       found2->second.end());
         _secondDerivMask[nParams * i + j] = canBeNonZero;
         _secondDerivMask[nParams * j + i] = canBeNonZero;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Report whether the second derivative with respect to parameters i and j
/// (indices in the space of all floatable parameters, matching Minuit's
/// external parameter indices) is identically zero because the parameters
/// share no additive term of the minimized function.
bool RooMinimizerFcn::secondDerivativeAlwaysVanishes(unsigned int i, unsigned int j) const
{
   std::call_once(_secondDerivMaskOnce, &RooMinimizerFcn::buildSecondDerivMask, this);
   return !_secondDerivMask[getNDim() * i + j];
}

/// \endcond
