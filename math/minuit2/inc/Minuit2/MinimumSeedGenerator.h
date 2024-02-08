// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MinimumSeedGenerator
#define ROOT_MyMinuit2_MinimumSeedGenerator

namespace ROOT {

namespace MyMinuit2 {

class MinimumSeed;
class MnFcn;
class GradientCalculator;
class MnUserParameterState;
class MnStrategy;
class AnalyticalGradientCalculator;

/** base class for seed generators (starting values); the seed generator
    prepares initial starting values from the input (MnUserParameterState)
    for the minimization;
 */

class MinimumSeedGenerator {

public:
   virtual ~MinimumSeedGenerator() {}

   virtual MinimumSeed
   operator()(const MnFcn &, const GradientCalculator &, const MnUserParameterState &, const MnStrategy &) const = 0;

   virtual MinimumSeed operator()(const MnFcn &, const AnalyticalGradientCalculator &, const MnUserParameterState &,
                                  const MnStrategy &) const = 0;
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_MinimumSeedGenerator
