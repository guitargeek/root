/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef ROOT_ROOFIT_TESTSTATISTICS_RooBinnedL
#define ROOT_ROOFIT_TESTSTATISTICS_RooBinnedL

#include <RooFit/TestStatistics/RooAbsL.h>
#include "RooAbsReal.h"

#include "Math/Util.h" // KahanSum

#include <vector>

// forward declarations
class RooAbsPdf;
class RooAbsData;
class RooChangeTracker;

namespace RooFit {
namespace TestStatistics {

class RooBinnedL : public RooAbsL {
public:
   RooBinnedL(RooAbsPdf *pdf, RooAbsData *data,
              RooFit::ZeroPredictionMode zeroPredMode = RooFit::ZeroPredictionMode::NaN, double zeroPredDelta = 1e-4);
   ~RooBinnedL() override;
   ROOT::Math::KahanSum<double>
   evaluatePartition(Section bins, std::size_t components_begin, std::size_t components_end) override;

   std::string GetClassName() const override { return "RooBinnedL"; }

private:
   mutable bool _first = true;        ///<!
   mutable std::vector<double> _binw; ///<!
   RooFit::ZeroPredictionMode const _zeroPredMode;
   double const _zeroPredDelta;
   mutable bool _zeroPredWarned = false; ///< whether we already warned about hitting the zero-prediction floor
   std::unique_ptr<RooChangeTracker> paramTracker_;
   Section lastSection_ = {0, 0}; // used for cache together with the parameter tracker
   mutable ROOT::Math::KahanSum<double> cachedResult_{0.};
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_RooBinnedL
