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

#include <RooFit/Detail/RooCategoryLookup.h>

#include <RooFit/EvalContext.h>

#include <algorithm>

/**
 * \class RooFit::Detail::RooCategoryLookup
 *
 * Piecewise-constant real-valued function of a category: evaluates to the
 * element of a list of real values selected by the current category state,
 * multiplied by a constant scale factor. The i-th list element corresponds to
 * the category state with index number i; states without a corresponding
 * element evaluate to zero.
 *
 * This is used to represent the penalty term that pdfs like RooMultiPdf add
 * to the negative log-likelihood (see RooAbsPdf::createCorrectionTerm()). The
 * penalty has to be a live function of the index category and not a constant,
 * so that it can influence the choice of the discrete parameters in the
 * discrete profiling loop of the RooMinimizer.
 */

namespace RooFit::Detail {

RooCategoryLookup::RooCategoryLookup(const char *name, const char *title, RooAbsCategory &cat, RooArgList const &values,
                                     double scale)
   : RooAbsReal(name, title),
     _cat("cat", "Input category", this, cat),
     _values("values", "Values, one for each category index number", this),
     _scale{scale}
{
   _values.add(values);
}

RooCategoryLookup::RooCategoryLookup(const RooCategoryLookup &other, const char *name)
   : RooAbsReal(other, name),
     _cat("cat", this, other._cat),
     _values("values", this, other._values),
     _scale{other._scale}
{
}

double RooCategoryLookup::evaluate() const
{
   const int idx = _cat->getCurrentIndex();
   if (idx < 0 || idx >= static_cast<int>(_values.size())) {
      return 0.0;
   }
   return _scale * static_cast<RooAbsReal const &>(_values[idx]).getVal();
}

void RooCategoryLookup::doEval(RooFit::EvalContext &ctx) const
{
   // The in-memory category state is authoritative: the evaluator itself
   // detects category changes via RooAbsCategory::getCurrentIndex().
   double val = 0.0;
   const int idx = _cat->getCurrentIndex();
   if (idx >= 0 && idx < static_cast<int>(_values.size())) {
      val = _scale * ctx.at(&_values[idx])[0];
   }
   std::span<double> output = ctx.output();
   std::fill(output.begin(), output.end(), val);
}

} // namespace RooFit::Detail
