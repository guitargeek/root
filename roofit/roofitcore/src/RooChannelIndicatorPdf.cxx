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

#include "RooFit/Detail/RooChannelIndicatorPdf.h"

#include "RooFit/EvalContext.h"

#include <TError.h>

#include <algorithm>

/**
 * \class RooFit::Detail::RooChannelIndicatorPdf
 *
 * Indicator density selecting one channel of a simultaneous fit, used by
 * RooSimultaneous::compileForNormSet() to represent the simultaneous pdf as
 * an ordinary mixture pdf. See the class documentation in the header for
 * details.
 */

namespace RooFit {
namespace Detail {

double RooChannelIndicatorPdf::evaluate() const
{
   for (std::size_t i = 0; i < _states.size(); ++i) {
      if (!matches(static_cast<RooAbsReal const &>(_indexVars[i]).getVal(), i)) {
         return 0.0;
      }
   }
   return 1.0;
}

void RooChannelIndicatorPdf::doEval(RooFit::EvalContext &ctx) const
{
   std::span<double> output = ctx.output();

   std::fill(output.begin(), output.end(), 1.0);
   for (std::size_t v = 0; v < _states.size(); ++v) {
      std::span<const double> indexVals = ctx.at(&_indexVars[v]);
      for (std::size_t i = 0; i < output.size(); ++i) {
         if (!matches(indexVals[indexVals.size() == 1 ? 0 : i], v)) {
            output[i] = 0.0;
         }
      }
   }
}

Int_t RooChannelIndicatorPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                                    const char * /*rangeName*/) const
{
   // Any subset of the index variables can be integrated analytically: the
   // code is a bitmask of the integrated variables (offset by one so that a
   // non-empty subset never maps to code zero). The range name is
   // deliberately ignored: the index stand-ins never carry named ranges in
   // the mixture compilation, because ranges on the index category drop the
   // excluded channels already at compile time. A named range that excluded
   // the target state would not be respected here.
   R__ASSERT(_indexVars.size() < 31);
   Int_t code = 0;
   for (std::size_t i = 0; i < _indexVars.size(); ++i) {
      if (allVars.find(_indexVars[i])) {
         analVars.add(_indexVars[i]);
         code |= (1 << i);
      }
   }
   return code == 0 ? 0 : code + 1;
}

double RooChannelIndicatorPdf::analyticalIntegral(Int_t code, const char * /*rangeName*/) const
{
   // Integrating an index variable with respect to the counting measure on
   // the channel index contributes an exact factor of one; the variables that
   // are not integrated contribute their indicator factor at the current
   // value.
   const Int_t mask = code - 1;
   double result = 1.0;
   for (std::size_t i = 0; i < _states.size(); ++i) {
      if (!(mask & (1 << i)) && !matches(static_cast<RooAbsReal const &>(_indexVars[i]).getVal(), i)) {
         result = 0.0;
      }
   }
   return result;
}

} // namespace Detail
} // namespace RooFit
