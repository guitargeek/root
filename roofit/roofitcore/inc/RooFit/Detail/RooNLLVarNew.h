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

#ifndef RooFit_RooNLLVarNew_h
#define RooFit_RooNLLVarNew_h

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooDataHist.h>
#include <RooGlobalFunc.h>
#include <RooListProxy.h>
#include <RooSetProxy.h>
#include <RooTemplateProxy.h>

#include <Math/Util.h>

#include <memory>
#include <vector>

namespace RooFit {
namespace Detail {

/// Template pdf built from the fit dataset, used for the bin-by-bin
/// likelihood offsetting (the Offset("bin") option of createNLL()). With a
/// mask, only the unmasked rows enter the template histogram, so that a
/// simultaneous pdf compiled into a mixture can build one template per
/// channel, normalized over the channel's own weight sum. With
/// `scaleByWeightSum`, the template values are scaled by sumOfWeights/e,
/// which absorbs the weight-sum-dependent parts of the per-channel offset
/// extended term of the channel-splitting path (see the implementation of
/// doEval()).
class RooOffsetPdf : public RooAbsPdf {
public:
   RooOffsetPdf(const char *name, const char *title, RooArgSet const &observables, RooAbsReal &weightVar,
                RooAbsReal *mask = nullptr, bool scaleByWeightSum = false)
      : RooAbsPdf(name, title),
        _observables("!observables", "List of observables", this),
        _weightVar{"!weightVar", "weightVar", this, weightVar, true, false},
        _scaleByWeightSum{scaleByWeightSum}
   {
      for (RooAbsArg *obs : observables) {
         _observables.add(*obs);
      }
      if (mask) {
         _mask = std::make_unique<RooTemplateProxy<RooAbsReal>>("!mask", "mask", this, *mask, true, false);
      }
   }
   RooOffsetPdf(const RooOffsetPdf &other, const char *name = nullptr)
      : RooAbsPdf(other, name),
        _observables("!servers", this, other._observables),
        _weightVar{"!weightVar", this, other._weightVar},
        _scaleByWeightSum{other._scaleByWeightSum}
   {
      if (other._mask) {
         _mask = std::make_unique<RooTemplateProxy<RooAbsReal>>("!mask", this, *other._mask);
      }
   }
   TObject *clone(const char *newname) const override { return new RooOffsetPdf(*this, newname); }

   /// Point the weight proxy to the weight variable of the likelihood; used
   /// when the offset pdf was created before the likelihood, by the
   /// simultaneous mixture compilation.
   void setWeightVar(RooAbsReal &weightVar) { _weightVar.setArg(weightVar); }

   void doEval(RooFit::EvalContext &ctx) const override;

private:
   double evaluate() const override { return 0.0; } // should never be called

   RooSetProxy _observables;
   RooTemplateProxy<RooAbsReal> _weightVar;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _mask;
   bool _scaleByWeightSum = false;
};

/// Per-row bin volumes for the concatenated binned likelihood of a
/// simultaneous mixture whose channels don't (all) use RooBinWidthFunctions
/// (e.g. models from before ROOT 6.26): the compiled values of such channels
/// are probability densities, and the likelihood multiplies each row by the
/// volume of the corresponding bin (see RooNLLVarNew::doEvalBinnedL()). The
/// rows of a channel are its bins in order, selected by the channel
/// indicator mask; rows of channels without a width vector (whose values are
/// already yields) get a unit volume.
class RooMixtureBinVolumes : public RooAbsReal {
public:
   RooMixtureBinVolumes(const char *name, const char *title, RooArgList const &indicators,
                        std::vector<std::vector<double>> binWidths)
      : RooAbsReal(name, title),
        _indicators{"!indicators", "channel indicators", this},
        _binWidths{std::move(binWidths)}
   {
      _indicators.add(indicators);
   }
   RooMixtureBinVolumes(const RooMixtureBinVolumes &other, const char *name = nullptr)
      : RooAbsReal(other, name), _indicators{"!indicators", this, other._indicators}, _binWidths{other._binWidths}
   {
   }
   TObject *clone(const char *newname) const override { return new RooMixtureBinVolumes(*this, newname); }

   void doEval(RooFit::EvalContext &ctx) const override;

private:
   double evaluate() const override { return _value; } // should never be called

   RooListProxy _indicators;
   std::vector<std::vector<double>> _binWidths; ///< per indicator, empty for unit volumes
};

/// Sum of the event weights of one channel of a simultaneous mixture,
/// selected by the channel indicator mask. Used as the per-channel expected
/// count in non-extended chi-squared fits of simultaneous pdfs compiled into
/// a mixture, where it plays the role that the dataset weight sum plays for
/// a single-channel chi-squared (see RooNLLVarNew::doEvalChi2()). The weight
/// variable is a placeholder until the chi-squared likelihood connects its
/// own weight variable via setWeightVar().
class RooChannelWeightSum : public RooAbsReal {
public:
   RooChannelWeightSum(const char *name, const char *title, RooAbsReal &mask, RooAbsReal &weightVar)
      : RooAbsReal(name, title),
        _mask{"!mask", "mask", this, mask, true, false},
        _weightVar{"!weightVar", "weightVar", this, weightVar, true, false}
   {
   }
   RooChannelWeightSum(const RooChannelWeightSum &other, const char *name = nullptr)
      : RooAbsReal(other, name), _mask{"!mask", this, other._mask}, _weightVar{"!weightVar", this, other._weightVar}
   {
   }
   TObject *clone(const char *newname) const override { return new RooChannelWeightSum(*this, newname); }

   /// Point the weight proxy to the weight variable of the likelihood.
   void setWeightVar(RooAbsReal &weightVar) { _weightVar.setArg(weightVar); }

   RooAbsReal const &mask() const { return *_mask; }
   RooAbsReal const &weightVar() const { return *_weightVar; }

   bool isReducerNode() const override { return true; }
   void doEval(RooFit::EvalContext &ctx) const override;

private:
   double evaluate() const override { return _value; } // should never be called

   RooTemplateProxy<RooAbsReal> _mask;
   RooTemplateProxy<RooAbsReal> _weightVar;

   ClassDefOverride(RooFit::Detail::RooChannelWeightSum, 0);
};

class RooNLLVarNew : public RooAbsReal {

public:
   // The names for the special variables that the RooNLLVarNew expects
   static constexpr const char *weightVarName = "_weight";
   static constexpr const char *weightVarNameSumW2 = "_weight_sumW2";
   static constexpr const char *binVolumeVarName = "_bin_volume";
   static constexpr const char *weightErrorLoVarName = "_weight_err_lo";
   static constexpr const char *weightErrorHiVarName = "_weight_err_hi";

   enum class Statistic {
      NLL,
      Chi2
   };

   /// Configuration struct for the unified constructor. Note that `offsetMode`
   /// only applies to `Statistic::NLL`, and `chi2ErrorType` only applies to
   /// `Statistic::Chi2`.
   struct Config {
      Statistic statistic = Statistic::NLL;
      bool extended = false;
      RooFit::OffsetMode offsetMode = RooFit::OffsetMode::None;
      RooDataHist::ErrorType chi2ErrorType = RooDataHist::Expected;
   };

   RooNLLVarNew(const char *name, const char *title, RooAbsReal &func, RooArgSet const &observables, Config const &cfg);
   RooNLLVarNew(const RooNLLVarNew &other, const char *name = nullptr);
   TObject *clone(const char *newname) const override { return new RooNLLVarNew(*this, newname); }

   /// Return default level for MINUIT error analysis.
   double defaultErrorLevel() const override { return _statistic == Statistic::Chi2 ? 1.0 : 0.5; }

   void doEval(RooFit::EvalContext &) const override;
   bool canComputeBatchWithCuda() const override { return _statistic == Statistic::NLL && !_binnedL && !_mixedBinnedL; }
   bool isReducerNode() const override { return true; }

   void applyWeightSquared(bool flag) override;

   void enableOffsetting(bool) override;

   void enableBinOffsetting(bool on = true) { _doBinOffset = on; }

   void setSimCount(int simCount) { _simCount = simCount; }

   enum class FuncMode {
      Pdf,
      ExtendedPdf,
      Function
   };

   RooAbsReal const &func() const { return *_func; }
   RooAbsReal const &weightVar() const { return *_weightVar; }
   RooAbsReal const &weightSquaredVar() const { return *_weightSquaredVar; }
   bool binnedL() const { return _binnedL; }
   bool mixedBinnedL() const { return _mixedBinnedL; }
   bool expectedEventsFolded() const { return _expectedEventsFolded; }
   int simCount() const { return _simCount; }
   Statistic statistic() const { return _statistic; }
   FuncMode funcMode() const { return _funcMode; }
   RooDataHist::ErrorType chi2ErrorType() const { return _chi2ErrorType; }
   RooAbsReal const *expectedEvents() const { return _expectedEvents ? &**_expectedEvents : nullptr; }
   RooAbsReal const *binnedRowsMask() const { return _binnedRowsMask ? &**_binnedRowsMask : nullptr; }
   RooAbsReal const *mixtureBinVolumes() const { return _mixtureBinVolumes ? &**_mixtureBinVolumes : nullptr; }
   RooAbsReal const *binVolumes() const { return _binVolumes ? &**_binVolumes : nullptr; }
   RooAbsReal const *weightErrLo() const { return _weightErrLo ? &**_weightErrLo : nullptr; }
   RooAbsReal const *weightErrHi() const { return _weightErrHi ? &**_weightErrHi : nullptr; }

private:
   double evaluate() const override { return _value; }
   double sumOfWeights(RooFit::EvalContext &, std::span<const double> weights, bool squared) const;
   void finalizeResult(RooFit::EvalContext &, ROOT::Math::KahanSum<double> result, double weightSum) const;
   void fillBinWidthsFromPdfBoundaries(RooAbsReal const &pdf, RooArgSet const &observables);
   void doEvalBinnedL(RooFit::EvalContext &, std::span<const double> preds, std::span<const double> weights) const;
   void doEvalMixed(RooFit::EvalContext &, std::span<const double> preds, std::span<const double> weights,
                    std::span<const double> weightsSumW2) const;
   void doEvalChi2(RooFit::EvalContext &, std::span<const double> preds, std::span<const double> weights,
                   std::span<const double> weightsSumW2) const;

   RooTemplateProxy<RooAbsReal> _func;
   RooTemplateProxy<RooAbsReal> _weightVar;
   RooTemplateProxy<RooAbsReal> _weightSquaredVar;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _expectedEvents;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _binnedRowsMask;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _mixtureBinVolumes;
   std::unique_ptr<RooTemplateProxy<RooAbsPdf>> _offsetPdf;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _binVolumes;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _weightErrLo;
   std::unique_ptr<RooTemplateProxy<RooAbsReal>> _weightErrHi;
   bool _weightSquared = false;
   bool _binnedL = false;
   bool _mixedBinnedL = false;
   /// Whether the expected-events proxy holds a gated-sum mixture total that
   /// is added to the likelihood directly (see the constructor).
   bool _expectedEventsFolded = false;
   bool _doOffset = false;
   bool _doBinOffset = false;
   Statistic _statistic = Statistic::NLL;
   FuncMode _funcMode = FuncMode::Pdf;
   RooDataHist::ErrorType _chi2ErrorType = RooDataHist::Expected;
   int _simCount = 1;
   std::vector<double> _binw;
   mutable ROOT::Math::KahanSum<double> _offset{0.}; ///<! Offset as KahanSum to avoid loss of precision
   mutable std::size_t _sumWeightGen = 0;            ///<! Input data generation of the cached weight sum
   mutable double _sumWeightCache = 0.0;             ///<! Cached sum of the event weights
   mutable std::size_t _sumWeight2Gen = 0;           ///<! Input data generation of the cached squared-weight sum
   mutable double _sumWeight2Cache = 0.0;            ///<! Cached sum of the squared event weights

   /// Cache for incremental evaluation, used only when the evaluator
   /// provides change tracking for the pdf values (see
   /// RooFit::EvalContext::changeTrackingEnabled()): the NLL reduction is
   /// done in fixed-size chunks whose partial results are cached, so that
   /// only the chunks in which the pdf values changed have to be reduced
   /// again. This makes the cost of a single-parameter variation in e.g. a
   /// simultaneous-fit mixture proportional to the size of the affected
   /// channel instead of the full dataset.
   struct ChunkCache {
      const double *probasPtr = nullptr;
      const double *weightsPtr = nullptr;
      std::size_t nEvents = 0;
      std::vector<double> sums;
      std::vector<double> carrys;
      std::vector<std::size_t> counts; ///< 3 entries per chunk: infinite, non-positive, NaN (1 per chunk for binned)
      std::vector<double> weightSums;  ///< per-chunk weight sums (binned likelihood only)
   };
   mutable ChunkCache _chunkCache; ///<!

   ClassDefOverride(RooFit::Detail::RooNLLVarNew, 0);
};

} // namespace Detail
} // namespace RooFit

#endif

/// \endcond
