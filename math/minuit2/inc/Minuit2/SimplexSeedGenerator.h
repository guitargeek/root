// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_SimplexSeedGenerator
#define ROOT_MyMinuit2_SimplexSeedGenerator

#include "MyMinuit2/MinimumSeedGenerator.h"

namespace ROOT {

namespace MyMinuit2 {

class MinimumSeed;
class MnFcn;
class MnUserParameterState;
class MnStrategy;

/**
   generate Simplex starting point (state)
 */
class SimplexSeedGenerator : public MinimumSeedGenerator {

public:
   SimplexSeedGenerator() {}

   ~SimplexSeedGenerator() override {}

   MinimumSeed
   operator()(const MnFcn &, const GradientCalculator &, const MnUserParameterState &, const MnStrategy &) const override;

   MinimumSeed operator()(const MnFcn &, const AnalyticalGradientCalculator &, const MnUserParameterState &,
                                  const MnStrategy &) const override;

private:
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_SimplexSeedGenerator
