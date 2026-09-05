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
#include <RooListProxy.h>

#include <TError.h>

#include <cmath>
#include <vector>

namespace RooFit {
namespace Detail {

/// Indicator density that selects one channel of a simultaneous fit.
///
/// Evaluates to one if every real-valued channel index variable is equal to
/// its given state index (within half a unit, since the index variables take
/// integer values loaded from category data columns), and to zero otherwise.
/// A simultaneous pdf with a plain index category has one index variable; one
/// whose index is a RooSuperCategory (e.g. from combining simultaneous pdfs)
/// has one index variable per input category.
///
/// The indicator is a probability density with respect to the *counting
/// measure* on the channel indices: summed over the possible index values it
/// gives exactly one. It therefore reports itself as self-normalized and
/// implements exact analytical integrals over any subset of the index
/// variables, so no numeric normalization object is ever attached to it.
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
      : RooChannelIndicatorPdf(name, title, RooArgList{indexVar}, std::vector<int>{state})
   {
   }

   RooChannelIndicatorPdf(const char *name, const char *title, RooArgList const &indexVars,
                          std::vector<int> const &states)
      : RooAbsPdf(name, title), _indexVars("indexVars", "channel index variables", this), _states{states}
   {
      R__ASSERT(indexVars.size() == states.size());
      _indexVars.add(indexVars);
   }

   RooChannelIndicatorPdf(const RooChannelIndicatorPdf &other, const char *name = nullptr)
      : RooAbsPdf(other, name), _indexVars("indexVars", this, other._indexVars), _states{other._states}
   {
   }

   TObject *clone(const char *newname) const override { return new RooChannelIndicatorPdf(*this, newname); }

   /// The density sums to one over the channel indices by construction.
   bool selfNormalized() const override { return true; }

   /// Numeric integration over the discontinuous indicator would be both
   /// wasteful and inexact; the analytical integrals below are exact.
   bool forceAnalyticalInt(const RooAbsArg & /*dep*/) const override { return true; }

   Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char *rangeName) const override;
   double analyticalIntegral(Int_t code, const char *rangeName) const override;

   std::vector<int> const &states() const { return _states; }
   RooArgList const &indexVars() const { return _indexVars; }

   void doEval(RooFit::EvalContext &) const override;

protected:
   double evaluate() const override;

private:
   bool matches(double index, std::size_t i) const { return std::abs(index - _states[i]) < 0.5; }

   RooListProxy _indexVars;
   std::vector<int> _states;

   ClassDefOverride(RooFit::Detail::RooChannelIndicatorPdf, 0);
};

} // namespace Detail
} // namespace RooFit

#endif
