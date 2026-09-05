/// \cond ROOFIT_INTERNAL

/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2021
 *   Emmanouil Michalainas, CERN 2021
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file RooNLLVarNew.cxx
\class RooNLLVarNew
\ingroup Roofitcore

This is a simple class designed to produce the nll values needed by the fitter.
This class calls functions from `RooBatchCompute` library to provide faster
computation times.
**/

#include "RooFit/Detail/RooNLLVarNew.h"

#include <RooHistPdf.h>
#include <RooBatchCompute.h>
#include <RooDataHist.h>
#include <RooNaNPacker.h>
#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooSetProxy.h>
#include <RooFit/Detail/MathFuncs.h>

#include "RooFitImplHelpers.h"

#include <ROOT/StringUtils.hxx>

#include <TClass.h>
#include <TMath.h>
#include <Math/Util.h>

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace RooFit::Detail {

// Declare constexpr static members to make them available if odr-used in C++14.
constexpr const char *RooNLLVarNew::weightVarName;
constexpr const char *RooNLLVarNew::weightVarNameSumW2;
constexpr const char *RooNLLVarNew::binVolumeVarName;
constexpr const char *RooNLLVarNew::weightErrorLoVarName;
constexpr const char *RooNLLVarNew::weightErrorHiVarName;

namespace {

// Use RooConstVar for dummies such that they don't get included in getParameters().
std::unique_ptr<RooConstVar> dummyVar(const char *name)
{
   return std::make_unique<RooConstVar>(name, name, 1.0);
}

} // namespace

void RooMixtureBinVolumes::doEval(RooFit::EvalContext &ctx) const
{
   std::span<double> output = ctx.output();

   // Rows that belong to no channel with a width vector keep a unit volume.
   std::fill(output.begin(), output.end(), 1.0);

   for (std::size_t c = 0; c < _binWidths.size(); ++c) {
      if (_binWidths[c].empty()) {
         continue;
      }
      std::span<const double> mask = ctx.at(static_cast<RooAbsReal const *>(&_indicators[c]));
      std::size_t bin = 0;
      for (std::size_t i = 0; i < output.size(); ++i) {
         if (mask[mask.size() == 1 ? 0 : i] > 0.5) {
            // The rows of a channel are its bins in order, like in the
            // per-channel likelihoods of the channel-splitting path.
            output[i] = _binWidths[c][std::min(bin, _binWidths[c].size() - 1)];
            ++bin;
         }
      }
   }
}

void RooChannelWeightSum::doEval(RooFit::EvalContext &ctx) const
{
   std::span<const double> mask = ctx.at(_mask);
   std::span<const double> weights = ctx.at(_weightVar);

   ROOT::Math::KahanSum<double> sum;
   const std::size_t n = std::max(mask.size(), weights.size());
   for (std::size_t i = 0; i < n; ++i) {
      if (mask[mask.size() == 1 ? 0 : i] > 0.5) {
         sum += weights[weights.size() == 1 ? 0 : i];
      }
   }
   ctx.output()[0] = sum.Sum();
}

void RooOffsetPdf::doEval(RooFit::EvalContext &ctx) const
{
   std::span<double> output = ctx.output();
   std::size_t nEvents = output.size();

   std::span<const double> weights = ctx.at(_weightVar);
   std::span<const double> mask = _mask ? ctx.at(*_mask) : std::span<const double>{};

   // Create the template histogram from the data. This operation is very
   // expensive, but since the offset only depends on the observables it
   // only has to be done once. Rows excluded by the mask (e.g. the rows of
   // foreign channels in a simultaneous mixture) don't enter the template.

   RooDataHist dataHist{"data", "data", _observables};
   // Loop over events to fill the histogram
   for (std::size_t i = 0; i < nEvents; ++i) {
      if (!mask.empty() && mask[mask.size() == 1 ? 0 : i] < 0.5) {
         continue;
      }
      for (auto *var : static_range_cast<RooRealVar *>(_observables)) {
         var->setVal(ctx.at(var)[i]);
      }
      dataHist.add(_observables, weights[weights.size() == 1 ? 0 : i]);
   }

   // Scaling the template by sumOfWeights/e reproduces, together with the
   // yield-scaled rows of an extended gated-sum mixture, exactly the
   // per-channel offset extended terms of the channel-splitting path,
   // (nu - W) - W*(log(nu) - log(W)) per channel: the row sum of the scaled
   // logarithms contributes sum_i w_i*log(T) + W*log(W) - W.
   const double scale = _scaleByWeightSum ? dataHist.sumEntries() / std::exp(1.0) : 1.0;

   // Lookup bin weights via RooHistPdf
   RooHistPdf pdf{"offsetPdf", "offsetPdf", _observables, dataHist};
   for (std::size_t i = 0; i < nEvents; ++i) {
      for (auto *var : static_range_cast<RooRealVar *>(_observables)) {
         var->setVal(ctx.at(var)[i]);
      }
      output[i] = scale * pdf.getVal(_observables);
   }
}

/// Construct either an NLL or a chi-squared test statistic.
/// \param func The pdf or function to evaluate. For `Statistic::NLL` a
/// RooAbsPdf is required, and for `Statistic::Chi2` any RooAbsReal is accepted.
RooNLLVarNew::RooNLLVarNew(const char *name, const char *title, RooAbsReal &func, RooArgSet const &observables,
                           Config const &cfg)
   : RooAbsReal(name, title),
     _func{"func", "func", this, func},
     _weightVar{"weightVar", "weightVar", this, dummyVar(weightVarName)},
     _weightSquaredVar{weightVarNameSumW2, weightVarNameSumW2, this, dummyVar("weightSquardVar")},
     _statistic{cfg.statistic},
     _chi2ErrorType{cfg.chi2ErrorType}
{
   auto *pdf = dynamic_cast<RooAbsPdf *>(&func);

   if (_statistic == Statistic::Chi2) {
      // Signal to RooEvaluatorWrapper::setData that zero-weight bins must be
      // retained (chi2 needs every bin's prediction, even where the data is
      // empty).
      setAttribute("Chi2EvaluationActive");
      _funcMode = !pdf ? FuncMode::Function : (cfg.extended ? FuncMode::ExtendedPdf : FuncMode::Pdf);
      if (pdf && pdf->getAttribute("MixtureChi2Active")) {
         // A simultaneous pdf compiled into a chi-squared mixture folds the
         // per-channel normalization factors (the expected channel yields,
         // or the per-channel data weight sums via RooChannelWeightSum
         // nodes) into its rows, so the values are used directly. The
         // weight-summing nodes still need this likelihood's weight
         // variable.
         _funcMode = FuncMode::Function;
         std::unique_ptr<RooArgSet> components{pdf->getComponents()};
         for (RooAbsArg *component : *components) {
            if (auto *weightSum = dynamic_cast<RooChannelWeightSum *>(component)) {
               weightSum->setWeightVar(*_weightVar);
            }
         }
      }
   } else {
      _binnedL = pdf && pdf->getAttribute("BinnedLikelihoodActive");
      _mixedBinnedL = pdf && pdf->getAttribute("MixedBinnedLikelihoodActive");
   }

   RooArgSet obs;
   func.getObservables(&observables, obs);

   const bool foldedExpectedEvents =
      _statistic == Statistic::NLL && pdf && pdf->getAttribute("MixtureFoldedExtendedEvents");

   if (_mixedBinnedL || foldedExpectedEvents) {
      // A pdf compiled from a simultaneous pdf into a gated-sum mixture
      // declares dedicated nodes in its computation graph: a mask selecting
      // the data rows of the binned-likelihood channels (when there are
      // both binned and unbinned channels), and (in extended fits) the total
      // expected events of the extendable channels. The latter is added to
      // the likelihood directly, because the rows of the gated sum already
      // carry the per-channel -log(expected yield) parts of the extended
      // terms. See RooSimultaneous::compileForNormSet().
      std::unique_ptr<RooArgSet> components{pdf->getComponents()};
      for (RooAbsArg *component : *components) {
         auto *asReal = dynamic_cast<RooAbsReal *>(component);
         if (!asReal) {
            continue;
         }
         if (component->getAttribute("MixtureBinnedRowsMask")) {
            _binnedRowsMask =
               std::make_unique<RooTemplateProxy<RooAbsReal>>("binnedRowsMask", "binnedRowsMask", this, *asReal);
         } else if (cfg.extended && component->getAttribute("MixtureExpectedEventsTotal")) {
            _expectedEvents =
               std::make_unique<RooTemplateProxy<RooAbsReal>>("expectedEvents", "expectedEvents", this, *asReal);
            _expectedEventsFolded = true;
         }
      }
      if (_mixedBinnedL && !_binnedRowsMask) {
         throw std::runtime_error("RooNLLVarNew: the mixed binned likelihood pdf declares no binned-rows mask");
      }
   }

   // Extended mode needs an expected-events function for both NLL and chi2
   // (NLL adds it as an extra additive term, chi2 uses it as the predicted
   // yield normalisation). Skip it for binned NLL (where the yields come
   // directly from the pdf), for the gated-sum mixtures (where the
   // expected-events sum node was already picked up above, or is
   // deliberately absent when no channel is extendable), and for chi2
   // Function mode (where the function values are directly the predicted
   // yields).
   const bool wantsExpectedEvents =
      pdf && ((_statistic == Statistic::NLL && cfg.extended && !_binnedL && !_mixedBinnedL && !foldedExpectedEvents) ||
              (_statistic == Statistic::Chi2 && _funcMode == FuncMode::ExtendedPdf));
   if (wantsExpectedEvents) {
      std::unique_ptr<RooAbsReal> expectedEvents = pdf->createExpectedEventsFunc(&obs);
      if (expectedEvents) {
         _expectedEvents =
            std::make_unique<RooTemplateProxy<RooAbsReal>>("expectedEvents", "expectedEvents", this, *expectedEvents);
         addOwnedComponents(std::move(expectedEvents));
      }
   }

   if (_statistic == Statistic::NLL) {
      // In the "BinnedLikelihoodActiveYields" mode, the pdf values can
      // directly be interpreted as yields and don't need to be multiplied by
      // the bin widths. That's why we don't need to even fill them in this
      // case. A simultaneous pdf compiled into a binned mixture provides its
      // per-row bin volumes as a node in its computation graph instead, when
      // not all of its channels are in yields mode.
      if (_binnedL) {
         std::unique_ptr<RooArgSet> components{pdf->getComponents()};
         for (RooAbsArg *component : *components) {
            if (component->getAttribute("MixtureBinVolumes")) {
               _mixtureBinVolumes = std::make_unique<RooTemplateProxy<RooAbsReal>>(
                  "mixtureBinVolumes", "mixtureBinVolumes", this, static_cast<RooAbsReal &>(*component));
            }
         }
         if (!_mixtureBinVolumes && !pdf->getAttribute("BinnedLikelihoodActiveYields")) {
            fillBinWidthsFromPdfBoundaries(*pdf, obs);
         }
      }

      enableOffsetting(cfg.offsetMode == RooFit::OffsetMode::Initial);
      enableBinOffsetting(cfg.offsetMode == RooFit::OffsetMode::Bin);

      // In the binned likelihood code path, we directly use that data weights
      // for the offsetting.
      if (!_binnedL && _doBinOffset) {
         // A simultaneous pdf compiled into a mixture provides its own
         // per-channel offset templates (see
         // RooSimultaneous::compileForNormSet()); discover them and connect
         // this likelihood's weight variable to them. Otherwise, build the
         // template pdf over all rows here.
         RooAbsPdf *mixtureOffset = nullptr;
         if (pdf) {
            std::unique_ptr<RooArgSet> components{pdf->getComponents()};
            for (RooAbsArg *component : *components) {
               if (component->getAttribute("MixtureBinOffsetPdf")) {
                  mixtureOffset = static_cast<RooAbsPdf *>(component);
               } else if (auto *channelOffset = dynamic_cast<RooOffsetPdf *>(component)) {
                  channelOffset->setWeightVar(*_weightVar);
               }
            }
         }
         if (mixtureOffset) {
            _offsetPdf = std::make_unique<RooTemplateProxy<RooAbsPdf>>("offsetPdf", "offsetPdf", this, *mixtureOffset);
         } else {
            auto offsetPdf = std::make_unique<RooOffsetPdf>("_offset_func", "_offset_func", obs, *_weightVar);
            _offsetPdf = std::make_unique<RooTemplateProxy<RooAbsPdf>>("offsetPdf", "offsetPdf", this, *offsetPdf);
            addOwnedComponents(std::move(offsetPdf));
         }
      }
   } else {
      // Chi2-only proxies: per-bin volumes, plus per-bin asymmetric errors
      // when Poisson error mode is requested.
      auto binVolumeDummy = std::make_unique<RooConstVar>(binVolumeVarName, binVolumeVarName, 1.0);
      _binVolumes = std::make_unique<RooTemplateProxy<RooAbsReal>>(binVolumeVarName, binVolumeVarName, this,
                                                                   *binVolumeDummy, true, false);
      addOwnedComponents(std::move(binVolumeDummy));

      if (_chi2ErrorType == RooDataHist::Poisson) {
         auto errLoDummy = std::make_unique<RooConstVar>(weightErrorLoVarName, weightErrorLoVarName, 1.0);
         auto errHiDummy = std::make_unique<RooConstVar>(weightErrorHiVarName, weightErrorHiVarName, 1.0);
         _weightErrLo = std::make_unique<RooTemplateProxy<RooAbsReal>>(weightErrorLoVarName, weightErrorLoVarName, this,
                                                                       *errLoDummy, true, false);
         _weightErrHi = std::make_unique<RooTemplateProxy<RooAbsReal>>(weightErrorHiVarName, weightErrorHiVarName, this,
                                                                       *errHiDummy, true, false);
         addOwnedComponents(std::move(errLoDummy));
         addOwnedComponents(std::move(errHiDummy));
      }
   }

   resetWeightVarNames();
}

RooNLLVarNew::RooNLLVarNew(const RooNLLVarNew &other, const char *name)
   : RooAbsReal(other, name),
     _func{"func", this, other._func},
     _weightVar{"weightVar", this, other._weightVar},
     _weightSquaredVar{"weightSquaredVar", this, other._weightSquaredVar},
     _weightSquared{other._weightSquared},
     _binnedL{other._binnedL},
     _mixedBinnedL{other._mixedBinnedL},
     _expectedEventsFolded{other._expectedEventsFolded},
     _doOffset{other._doOffset},
     _doBinOffset{other._doBinOffset},
     _statistic{other._statistic},
     _funcMode{other._funcMode},
     _chi2ErrorType{other._chi2ErrorType},
     _simCount{other._simCount},
     _prefix{other._prefix},
     _binw{other._binw}
{
   if (other._expectedEvents) {
      _expectedEvents = std::make_unique<RooTemplateProxy<RooAbsReal>>("expectedEvents", this, *other._expectedEvents);
   }
   if (other._binnedRowsMask) {
      _binnedRowsMask = std::make_unique<RooTemplateProxy<RooAbsReal>>("binnedRowsMask", this, *other._binnedRowsMask);
   }
   if (other._mixtureBinVolumes) {
      _mixtureBinVolumes =
         std::make_unique<RooTemplateProxy<RooAbsReal>>("mixtureBinVolumes", this, *other._mixtureBinVolumes);
   }
   if (other._binVolumes) {
      _binVolumes = std::make_unique<RooTemplateProxy<RooAbsReal>>(binVolumeVarName, this, *other._binVolumes);
   }
   if (other._weightErrLo) {
      _weightErrLo = std::make_unique<RooTemplateProxy<RooAbsReal>>(weightErrorLoVarName, this, *other._weightErrLo);
   }
   if (other._weightErrHi) {
      _weightErrHi = std::make_unique<RooTemplateProxy<RooAbsReal>>(weightErrorHiVarName, this, *other._weightErrHi);
   }
}

void RooNLLVarNew::fillBinWidthsFromPdfBoundaries(RooAbsReal const &pdf, RooArgSet const &observables)
{
   // Check if the bin widths were already filled
   if (!_binw.empty()) {
      return;
   }

   if (observables.size() != 1) {
      throw std::runtime_error("BinnedPdf optimization only works with a 1D pdf.");
   } else {
      auto *var = static_cast<RooRealVar *>(observables.first());
      std::list<double> *boundaries = pdf.binBoundaries(*var, var->getMin(), var->getMax());
      std::list<double>::iterator biter = boundaries->begin();
      _binw.resize(boundaries->size() - 1);
      double lastBound = (*biter);
      ++biter;
      int ibin = 0;
      while (biter != boundaries->end()) {
         _binw[ibin] = (*biter) - lastBound;
         lastBound = (*biter);
         ibin++;
         ++biter;
      }
   }
}

void RooNLLVarNew::doEvalBinnedL(RooFit::EvalContext &ctx, std::span<const double> preds,
                                 std::span<const double> weights) const
{
   const bool predsAreYields = _binw.empty() && !_mixtureBinVolumes;
   std::span<const double> binVolumes = _mixtureBinVolumes ? ctx.at(*_mixtureBinVolumes) : std::span<const double>{};

   // If the evaluator tracks in which event range the pdf values changed,
   // reduce the binned likelihood in fixed-size chunks and cache the partial
   // results, so that only the chunks with changed yields have to be reduced
   // again (see the analogous code for the unbinned reduction in doEval()).
   // This makes the cost of a single-parameter variation in a concatenated
   // binned likelihood proportional to the size of the affected channel. The
   // chunks are much smaller than for the unbinned reduction, because binned
   // datasets have far fewer entries.
   if (ctx.changeTrackingEnabled() && preds.size() > 1 && weights.size() == preds.size()) {
      constexpr std::size_t chunkSize = 64;
      const std::size_t n = preds.size();
      const std::size_t nChunks = (n + chunkSize - 1) / chunkSize;
      ChunkCache &cache = _chunkCache;

      auto [changedBegin, changedEnd] = ctx.changedRange(&*_func);

      const bool rebuild = cache.probasPtr != preds.data() || cache.weightsPtr != weights.data() ||
                           cache.nEvents != n || cache.sums.size() != nChunks;
      if (rebuild) {
         cache.probasPtr = preds.data();
         cache.weightsPtr = weights.data();
         cache.nEvents = n;
         cache.sums.assign(nChunks, 0.0);
         cache.carrys.assign(nChunks, 0.0);
         cache.counts.assign(nChunks, 0);
         cache.weightSums.assign(nChunks, 0.0);
         changedBegin = 0;
         changedEnd = n;
      }

      const std::size_t firstChunk = changedBegin / chunkSize;
      const std::size_t endChunk = changedEnd == 0 ? 0 : (std::min(changedEnd, n) - 1) / chunkSize + 1;
      for (std::size_t c = firstChunk; c < endChunk; ++c) {
         const std::size_t begin = c * chunkSize;
         const std::size_t end = std::min(begin + chunkSize, n);
         ROOT::Math::KahanSum<double> chunkSum;
         ROOT::Math::KahanSum<double> chunkWeightSum;
         std::size_t nErrors = 0;
         for (std::size_t i = begin; i < end; ++i) {
            const double N = weights[i];
            double mu = preds[i];
            if (!predsAreYields) {
               mu *= binVolumes.empty() ? _binw[i] : binVolumes[i];
            }
            if (mu <= 0 && N > 0) {
               ++nErrors;
            } else {
               chunkSum += RooFit::Detail::MathFuncs::nll(mu, N, true, _doBinOffset);
               chunkWeightSum += N;
            }
         }
         cache.sums[c] = chunkSum.Sum();
         cache.carrys[c] = chunkSum.Carry();
         cache.counts[c] = nErrors;
         cache.weightSums[c] = chunkWeightSum.Sum();
      }

      ROOT::Math::KahanSum<double> total;
      double sumWeight = 0.0;
      std::size_t nErrors = 0;
      for (std::size_t c = 0; c < nChunks; ++c) {
         total += ROOT::Math::KahanSum<double>{cache.sums[c], cache.carrys[c]};
         sumWeight += cache.weightSums[c];
         nErrors += cache.counts[c];
      }
      // The error condition (data present where zero events are predicted)
      // has to be reported on every evaluation, also for the bins whose
      // cached chunks were not recomputed.
      for (std::size_t i = 0; i < nErrors; ++i) {
         logEvalError("Observed events in a bin with zero event yield");
      }

      finalizeResult(ctx, total, sumWeight);
      return;
   }

   ROOT::Math::KahanSum<double> result{0.0};
   ROOT::Math::KahanSum<double> sumWeightKahanSum{0.0};

   for (std::size_t i = 0; i < preds.size(); ++i) {

      // Calculate log(Poisson(N|mu) for this bin
      double N = weights[i];
      double mu = preds[i];
      if (!predsAreYields) {
         mu *= binVolumes.empty() ? _binw[i] : binVolumes[i];
      }

      if (mu <= 0 && N > 0) {
         // Catch error condition: data present where zero events are predicted
         logEvalError(Form("Observed %f events in bin %lu with zero event yield", N, (unsigned long)i));
      } else {
         result += RooFit::Detail::MathFuncs::nll(mu, N, true, _doBinOffset);
         sumWeightKahanSum += N;
      }
   }

   finalizeResult(ctx, result, sumWeightKahanSum.Sum());
}

/// Return the sum of the event weights (or of the squared event weights).
/// The sum only changes when new input data is loaded into the evaluation
/// context, which is tracked with EvalContext::inputGeneration(). Caching it
/// avoids a reduction for every evaluation, which in CUDA mode would also
/// synchronize the stream.
double RooNLLVarNew::sumOfWeights(RooFit::EvalContext &ctx, std::span<const double> weights, bool squared) const
{
   std::size_t &generation = squared ? _sumWeight2Gen : _sumWeightGen;
   double &cache = squared ? _sumWeight2Cache : _sumWeightCache;
   if (generation != ctx.inputGeneration()) {
      cache = RooBatchCompute::reduceSum(ctx.config(this), weights.data(), weights.size());
      generation = ctx.inputGeneration();
   }
   return cache;
}

/// Reduce the likelihood of a pdf compiled from a simultaneous pdf with both
/// binned-likelihood and unbinned channels. The rows of the concatenated
/// dataset are heterogeneous: rows of the binned channels are histogram bins
/// (with the observed counts as weights and the compiled pdf values being the
/// expected yields), rows of the unbinned channels are events (with the
/// compiled pdf values being probability densities, scaled by the expected
/// channel yields in extended fits). The mask provided by the compiled pdf
/// selects the binned rows, which contribute Poisson terms exactly like in
/// doEvalBinnedL(); the unbinned rows contribute -weight * log(value) with
/// the same conventions as RooBatchCompute::reduceNLL().
void RooNLLVarNew::doEvalMixed(RooFit::EvalContext &ctx, std::span<const double> preds, std::span<const double> weights,
                               std::span<const double> weightsSumW2) const
{
   std::span<const double> mask = ctx.at(*_binnedRowsMask);
   std::span<const double> weightSpan = _weightSquared ? weightsSumW2 : weights;

   ROOT::Math::KahanSum<double> result;
   ROOT::Math::KahanSum<double> sumWeightK;
   // Weight sums over the unbinned rows only, for the sum-of-weights-squared
   // scaling of the expected-events term in extended weighted fits.
   ROOT::Math::KahanSum<double> sumWeightUnbinned;
   ROOT::Math::KahanSum<double> sumWeight2Unbinned;
   std::size_t nBinnedErrors = 0;
   std::size_t nNonPositive = 0;
   std::size_t nInfinite = 0;
   std::size_t nNaN = 0;
   double badness = 0.0;

   for (std::size_t i = 0; i < preds.size(); ++i) {
      const double w = weightSpan[weightSpan.size() == 1 ? 0 : i];
      if (mask[i] > 0.5) {
         // Binned row: log(Poisson(N|mu)) with the pdf value as the yield.
         // Like in doEvalBinnedL(), bins with data in a zero-yield prediction
         // are excluded from the sum and from the weight sum, and reported as
         // evaluation errors.
         const double N = w;
         const double mu = preds[i];
         if (mu <= 0 && N > 0) {
            ++nBinnedErrors;
         } else {
            result += RooFit::Detail::MathFuncs::nll(mu, N, true, false);
            sumWeightK += N;
         }
         continue;
      }

      // Unbinned row: -weight * log(value), with the error handling of
      // RooBatchCompute::reduceNLL(), including the skipping of zero-weight
      // entries.
      if (w == 0.0) {
         continue;
      }
      const double p = preds[i];
      double term = 0.0;
      if (p <= 0.0) {
         ++nNonPositive;
         term = std::log(p);
         badness += -p;
      } else if (std::isnan(p)) {
         ++nNaN;
         term = p;
         badness += RooNaNPacker::unpackNaN(p);
      } else {
         if (std::isinf(p)) {
            ++nInfinite;
         }
         term = std::log(p);
      }
      result += -w * term;
      sumWeightK += w;
      sumWeightUnbinned += weights[weights.size() == 1 ? 0 : i];
      sumWeight2Unbinned += weightsSumW2[weightsSumW2.size() == 1 ? 0 : i];
   }

   for (std::size_t i = 0; i < nBinnedErrors; ++i) {
      logEvalError("Observed events in a bin with zero event yield");
   }
   if (nInfinite > 0) {
      oocoutW(&*_func, Eval) << "RooAbsPdf::getLogVal(" << _func->GetName()
                             << ") WARNING: top-level pdf has some infinite values" << std::endl;
   }
   for (std::size_t i = 0; i < nNonPositive; ++i) {
      _func->logEvalError("getLogVal() top-level p.d.f not greater than zero");
   }
   for (std::size_t i = 0; i < nNaN; ++i) {
      _func->logEvalError("getLogVal() top-level p.d.f evaluates to NaN");
   }
   if (badness != 0.0) {
      // Some events with evaluation errors: return the "badness" of the
      // errors, packed like in RooBatchCompute::reduceNLL().
      result = ROOT::Math::KahanSum<double>{RooNaNPacker::packFloatIntoNaN(badness), 0.0};
   }

   if (_expectedEvents) {
      // The unbinned rows are the negative logarithms of the yield-scaled
      // channel densities, so their sum already includes the
      // -sumWeight * log(expectedEvents) part of the extended term; only the
      // sum of the expected events of the unbinned channels remains. In
      // weighted fits with the weight-squared correction, the term is scaled
      // by the ratio of the weight sums like in RooAbsPdf::extendedTerm().
      double term = ctx.at(*_expectedEvents)[0];
      if (_weightSquared && sumWeightUnbinned.Sum() > 0.0) {
         term *= sumWeight2Unbinned.Sum() / sumWeightUnbinned.Sum();
      }
      result += term;
   }

   finalizeResult(ctx, result, sumWeightK.Sum());
}

void RooNLLVarNew::doEvalChi2(RooFit::EvalContext &ctx, std::span<const double> preds, std::span<const double> weights,
                              std::span<const double> weightsSumW2) const
{
   // Error type None implies zero sigma for every bin: the chi2 is undefined
   // everywhere but empty bins. Match the legacy behaviour of returning zero.
   if (_chi2ErrorType == RooDataHist::None) {
      finalizeResult(ctx, ROOT::Math::KahanSum<double>{0.0}, 0.0);
      return;
   }

   std::span<const double> binVol = ctx.at(*_binVolumes);
   std::span<const double> errLo = _weightErrLo ? ctx.at(*_weightErrLo) : std::span<const double>{};
   std::span<const double> errHi = _weightErrHi ? ctx.at(*_weightErrHi) : std::span<const double>{};

   const double sumWeight = sumOfWeights(ctx, weights, false);

   double normFactor = 1.0;
   switch (_funcMode) {
   case FuncMode::Pdf: normFactor = sumWeight; break;
   case FuncMode::ExtendedPdf: normFactor = ctx.at(*_expectedEvents)[0]; break;
   case FuncMode::Function: normFactor = 1.0; break;
   }

   ROOT::Math::KahanSum<double> result{0.0};
   ROOT::Math::KahanSum<double> sumWeightKahanSum{0.0};

   for (std::size_t i = 0; i < preds.size(); ++i) {
      const double N = weights[i];
      const double mu = preds[i] * normFactor * binVol[i];
      const double diff = mu - N;

      double sigma2;
      switch (_chi2ErrorType) {
      case RooDataHist::SumW2: sigma2 = weightsSumW2[i]; break;
      case RooDataHist::Poisson: {
         // Poisson errors are asymmetric: choose the side facing the prediction.
         const double err = diff > 0 ? errHi[i] : errLo[i];
         sigma2 = err * err;
         break;
      }
      default: sigma2 = mu; break; // Expected
      }

      // Skip bins where data, prediction and error are all zero (matches legacy RooChi2Var).
      if (sigma2 == 0.0 && N == 0.0 && mu == 0.0) {
         continue;
      }
      if (sigma2 <= 0.0) {
         logEvalError(Form("chi2 bin %lu has non-positive error; term replaced with NaN", (unsigned long)i));
         result += std::numeric_limits<double>::quiet_NaN();
         continue;
      }

      result += diff * diff / sigma2;
      sumWeightKahanSum += N;
   }

   finalizeResult(ctx, result, sumWeightKahanSum.Sum());
}

void RooNLLVarNew::doEval(RooFit::EvalContext &ctx) const
{
   std::span<const double> weights = ctx.at(_weightVar);
   std::span<const double> weightsSumW2 = ctx.at(_weightSquaredVar);

   if (_statistic == Statistic::Chi2) {
      return doEvalChi2(ctx, ctx.at(&*_func), weights, weightsSumW2);
   }

   if (_binnedL) {
      return doEvalBinnedL(ctx, ctx.at(&*_func), _weightSquared ? weightsSumW2 : weights);
   }

   if (_mixedBinnedL) {
      return doEvalMixed(ctx, ctx.at(&*_func), weights, weightsSumW2);
   }

   auto config = ctx.config(this);

   auto probas = ctx.at(_func);
   std::span<const double> weightSpan = _weightSquared ? weightsSumW2 : weights;

   // The weight sums only depend on the input data, so they are cached in the
   // RooNLLVarNew keyed on the evaluation context's input generation.
   double sumWeight = sumOfWeights(ctx, weights, false);
   double sumWeight2 = 0.;
   if (_expectedEvents && _weightSquared) {
      sumWeight2 = sumOfWeights(ctx, weightsSumW2, true);
   }

   // If the evaluator tracks in which event range the pdf values changed,
   // reduce the NLL in fixed-size chunks and cache the partial results, so
   // that only the chunks with changed pdf values have to be reduced again.
   const bool incremental = !config.useCuda() && !_doBinOffset && ctx.changeTrackingEnabled() && probas.size() > 1;

   RooBatchCompute::ReduceNLLOutput nllOut;

   if (!incremental) {
      nllOut = RooBatchCompute::reduceNLL(config, probas, weightSpan,
                                          _doBinOffset ? ctx.at(*_offsetPdf) : std::span<const double>{});
   } else {
      constexpr std::size_t chunkSize = 4096;
      const std::size_t n = probas.size();
      const std::size_t nChunks = (n + chunkSize - 1) / chunkSize;
      ChunkCache &cache = _chunkCache;

      auto [changedBegin, changedEnd] = ctx.changedRange(&*_func);

      const bool rebuild = cache.probasPtr != probas.data() || cache.weightsPtr != weightSpan.data() ||
                           cache.nEvents != n || cache.sums.size() != nChunks;
      if (rebuild) {
         cache.probasPtr = probas.data();
         cache.weightsPtr = weightSpan.data();
         cache.nEvents = n;
         cache.sums.assign(nChunks, 0.0);
         cache.carrys.assign(nChunks, 0.0);
         cache.counts.assign(3 * nChunks, 0);
         changedBegin = 0;
         changedEnd = n;
      }

      const std::size_t firstChunk = changedBegin / chunkSize;
      const std::size_t endChunk = changedEnd == 0 ? 0 : (std::min(changedEnd, n) - 1) / chunkSize + 1;
      for (std::size_t c = firstChunk; c < endChunk; ++c) {
         const std::size_t begin = c * chunkSize;
         const std::size_t len = std::min(chunkSize, n - begin);
         std::span<const double> probasChunk{probas.data() + begin, len};
         std::span<const double> weightsChunk =
            weightSpan.size() == 1 ? weightSpan : std::span<const double>{weightSpan.data() + begin, len};
         auto out = RooBatchCompute::reduceNLL(config, probasChunk, weightsChunk, {});
         cache.sums[c] = out.nllSum;
         cache.carrys[c] = out.nllSumCarry;
         cache.counts[3 * c] = out.nInfiniteValues;
         cache.counts[3 * c + 1] = out.nNonPositiveValues;
         cache.counts[3 * c + 2] = out.nNaNValues;
      }

      ROOT::Math::KahanSum<double> total;
      for (std::size_t c = 0; c < nChunks; ++c) {
         total += ROOT::Math::KahanSum<double>{cache.sums[c], cache.carrys[c]};
         nllOut.nInfiniteValues += cache.counts[3 * c];
         nllOut.nNonPositiveValues += cache.counts[3 * c + 1];
         nllOut.nNaNValues += cache.counts[3 * c + 2];
      }
      nllOut.nllSum = total.Sum();
      nllOut.nllSumCarry = total.Carry();
   }

   if (nllOut.nInfiniteValues > 0) {
      oocoutW(&*_func, Eval) << "RooAbsPdf::getLogVal(" << _func->GetName()
                             << ") WARNING: top-level pdf has some infinite values" << std::endl;
   }
   for (std::size_t i = 0; i < nllOut.nNonPositiveValues; ++i) {
      _func->logEvalError("getLogVal() top-level p.d.f not greater than zero");
   }
   for (std::size_t i = 0; i < nllOut.nNaNValues; ++i) {
      _func->logEvalError("getLogVal() top-level p.d.f evaluates to NaN");
   }

   if (_expectedEvents) {
      std::span<const double> expected = ctx.at(*_expectedEvents);
      if (_expectedEventsFolded) {
         // The rows of a gated-sum mixture already carry the per-channel
         // -weight * log(expected yield) parts of the extended terms, so
         // only the summed expected events remain. In weighted fits with
         // the weight-squared correction, the term is scaled by the ratio
         // of the weight sums like in RooAbsPdf::extendedTerm().
         double term = expected[0];
         if (_weightSquared && sumWeight > 0.0) {
            term *= sumWeight2 / sumWeight;
         }
         nllOut.nllSum += term;
      } else {
         // The unbinned NLL path is only reached for pdf inputs, so the cast is safe.
         auto &pdf = static_cast<RooAbsPdf &>(const_cast<RooAbsReal &>(*_func));
         nllOut.nllSum += pdf.extendedTerm(sumWeight, expected[0], _weightSquared ? sumWeight2 : 0.0, _doBinOffset);
      }
   }

   finalizeResult(ctx, {nllOut.nllSum, nllOut.nllSumCarry}, sumWeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the prefix for the special variables of this NLL, like weights or bin
/// volumes.
/// \param[in] prefix The prefix to add to the observables and weight names.
void RooNLLVarNew::setPrefix(std::string const &prefix)
{
   _prefix = prefix;

   resetWeightVarNames();
}

void RooNLLVarNew::resetWeightVarNames()
{
   _weightVar->SetName((_prefix + weightVarName).c_str());
   _weightSquaredVar->SetName((_prefix + weightVarNameSumW2).c_str());
   if (_offsetPdf && !(*_offsetPdf)->getAttribute("MixtureBinOffsetPdf")) {
      // Only the template pdf built by this class is renamed; a discovered
      // mixture offset node keeps its name in the compiled graph.
      (*_offsetPdf)->SetName((_prefix + "_offset_func").c_str());
   }
   if (_binVolumes) {
      (*_binVolumes)->SetName((_prefix + binVolumeVarName).c_str());
   }
   if (_weightErrLo) {
      (*_weightErrLo)->SetName((_prefix + weightErrorLoVarName).c_str());
   }
   if (_weightErrHi) {
      (*_weightErrHi)->SetName((_prefix + weightErrorHiVarName).c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Toggles the weight square correction.
void RooNLLVarNew::applyWeightSquared(bool flag)
{
   if (_statistic == Statistic::Chi2) {
      if (flag) {
         coutW(Fitting) << "RooNLLVarNew::applyWeightSquared(" << GetName()
                        << ") has no effect on a chi-squared evaluator; ignoring." << std::endl;
      }
      return;
   }
   _weightSquared = flag;
}

void RooNLLVarNew::enableOffsetting(bool flag)
{
   _doOffset = flag;
   _offset = ROOT::Math::KahanSum<double>{};
}

void RooNLLVarNew::finalizeResult(RooFit::EvalContext &ctx, ROOT::Math::KahanSum<double> result, double weightSum) const
{
   // If part of simultaneous PDF normalize probability over
   // number of simultaneous PDFs: -sum(log(p/n)) = -sum(log(p)) + N*log(n)
   // If we do bin-by bin offsetting, we don't do this because it cancels out.
   // The correction is specific to NLL; it has no meaning for chi2.
   if (_statistic == Statistic::NLL && !_doBinOffset && _simCount > 1) {
      result += weightSum * std::log(static_cast<double>(_simCount));
   }

   // Check if value offset flag is set.
   if (_doOffset) {

      // If no offset is stored enable this feature now
      if (_offset.Sum() == 0 && _offset.Carry() == 0 && (result.Sum() != 0 || result.Carry() != 0)) {
         _offset = result;
      }
   }
   ctx.setOutputWithOffset(this, result, _offset);
}

} // namespace RooFit::Detail

/// \endcond
