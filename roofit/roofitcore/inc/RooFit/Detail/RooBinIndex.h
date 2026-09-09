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

#ifndef RooFit_Detail_RooBinIndex_h
#define RooFit_Detail_RooBinIndex_h

#include <RooAbsBinning.h>
#include <RooAbsReal.h>
#include <RooListProxy.h>

#include <memory>
#include <vector>

class RooDataHist;

namespace RooFit::Detail {

class CompileContext;

/// A node in the computation graph that represents the flattened bin index
/// corresponding to the current values of a set of observables in a given
/// binned structure.
///
/// This class is an implementation detail of the computation graphs created
/// by RooAbsArg::compileForNormSet(). It is used to deduplicate bin index
/// calculations between histogram-based objects that share the same
/// observables and binnings, like the many template histograms in a
/// HistFactory channel. It is not meant to be used in public user-facing
/// models, and it doesn't support ROOT IO.
///
/// The represented index is `sum_i coefs[i] * binnings[i]->binNumber(vars[i])`.
/// Since the RooFit batch evaluation interface only supports real-valued
/// outputs, the bin index is represented as a double (integers are exactly
/// representable well beyond any realistic number of bins).
class RooBinIndex : public RooAbsReal {
public:
   RooBinIndex(const char *name, RooArgList const &vars, std::vector<RooAbsBinning const *> const &binnings,
               std::vector<int> const &coefs);
   RooBinIndex(const RooBinIndex &other, const char *name = nullptr);
   TObject *clone(const char *newname = nullptr) const override { return new RooBinIndex(*this, newname); }

   /// Get a bin index node for the given observables and binnings from the
   /// compilation context, creating and registering a new one if no
   /// equivalent node exists yet. The first creator passes ownership of the
   /// node to `owner`.
   static RooBinIndex *getOrCreate(CompileContext &ctx, RooAbsArg &owner, RooArgList const &vars,
                                   std::vector<RooAbsBinning const *> const &binnings, std::vector<int> const &coefs);

   /// Like getOrCreate(), for the binnings of a RooDataHist with the
   /// RooDataHist index convention (last variable is the fastest-running
   /// index). Returns `nullptr` if the histogram structure is not supported
   /// (e.g. category dimensions or parameterized binnings).
   static RooBinIndex *getOrCreateForDataHist(CompileContext &ctx, RooAbsArg &owner, RooArgList const &vars,
                                              RooDataHist const &dataHist);

   RooArgList const &vars() const { return _vars; }
   std::vector<std::unique_ptr<RooAbsBinning>> const &binnings() const { return _binnings; }
   std::vector<int> const &coefs() const { return _coefs; }

protected:
   double evaluate() const override;
   void doEval(RooFit::EvalContext &) const override;

private:
   RooListProxy _vars;
   std::vector<std::unique_ptr<RooAbsBinning>> _binnings;
   std::vector<int> _coefs;
   mutable std::vector<int> _intBuffer; ///<! working buffer for batch evaluations

   ClassDefOverride(RooFit::Detail::RooBinIndex, 0);
};

} // namespace RooFit::Detail

#endif
