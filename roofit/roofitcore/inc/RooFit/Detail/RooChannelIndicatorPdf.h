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

#ifndef RooFit_Detail_RooChannelIndicatorPdf_h
#define RooFit_Detail_RooChannelIndicatorPdf_h

#include <RooAbsPdf.h>
#include <RooRealProxy.h>

#include <cmath>

namespace RooFit {
namespace Detail {

/// Indicator density that selects one channel of a simultaneous fit.
///
/// Evaluates to one if the real-valued channel index variable is equal to the
/// given state index (within half a unit, since the index variable takes
/// integer values loaded from a category data column), and to zero otherwise.
///
/// The indicator is a probability density with respect to the *counting
/// measure* on the channel index: summed over the possible index values it
/// gives exactly one. It therefore reports itself as self-normalized and
/// implements an exact unit analytical integral over the index variable, so
/// no numeric normalization object is ever attached to it.
///
/// This class is used by RooSimultaneous::compileForNormSet() to represent
/// the simultaneous pdf as an ordinary mixture,
/// \f[
///   P(x, c) = \sum_s w_s \; \mathbf{1}[c = s] \; f_s(x),
/// \f]
/// so that the likelihood machinery downstream doesn't need any special
/// treatment of the simultaneous case.
class RooChannelIndicatorPdf : public RooAbsPdf {
public:
   RooChannelIndicatorPdf(const char *name, const char *title, RooAbsReal &indexVar, int state)
      : RooAbsPdf(name, title), _indexVar("indexVar", "channel index variable", this, indexVar), _state{state}
   {
   }

   RooChannelIndicatorPdf(const RooChannelIndicatorPdf &other, const char *name = nullptr)
      : RooAbsPdf(other, name), _indexVar("indexVar", this, other._indexVar), _state{other._state}
   {
   }

   TObject *clone(const char *newname) const override { return new RooChannelIndicatorPdf(*this, newname); }

   /// The density sums to one over the channel index by construction.
   bool selfNormalized() const override { return true; }

   /// Numeric integration over the discontinuous indicator would be both
   /// wasteful and inexact; the analytical integral below is exact.
   bool forceAnalyticalInt(const RooAbsArg & /*dep*/) const override { return true; }

   Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char *rangeName) const override;
   double analyticalIntegral(Int_t code, const char *rangeName) const override;

   int state() const { return _state; }
   RooAbsReal const &indexVar() const { return *_indexVar; }

   void doEval(RooFit::EvalContext &) const override;

protected:
   double evaluate() const override { return matches(_indexVar) ? 1.0 : 0.0; }

private:
   bool matches(double index) const { return std::abs(index - _state) < 0.5; }

   RooRealProxy _indexVar;
   int _state;

   ClassDefOverride(RooFit::Detail::RooChannelIndicatorPdf, 0);
};

} // namespace Detail
} // namespace RooFit

#endif
