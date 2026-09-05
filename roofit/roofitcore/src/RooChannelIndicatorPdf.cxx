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

void RooChannelIndicatorPdf::doEval(RooFit::EvalContext &ctx) const
{
   std::span<const double> indexVals = ctx.at(_indexVar);
   std::span<double> output = ctx.output();

   for (std::size_t i = 0; i < output.size(); ++i) {
      output[i] = matches(indexVals[indexVals.size() == 1 ? 0 : i]) ? 1.0 : 0.0;
   }
}

Int_t RooChannelIndicatorPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                                    const char * /*rangeName*/) const
{
   if (matchArgs(allVars, analVars, _indexVar)) {
      return 1;
   }
   return 0;
}

double RooChannelIndicatorPdf::analyticalIntegral(Int_t code, const char * /*rangeName*/) const
{
   R__ASSERT(code == 1);
   // Unit integral with respect to the counting measure on the channel index,
   // exact by construction.
   return 1.0;
}

} // namespace Detail
} // namespace RooFit
