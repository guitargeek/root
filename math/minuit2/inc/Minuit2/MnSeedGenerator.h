// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MnSeedGenerator
#define ROOT_MyMinuit2_MnSeedGenerator

#include "MyMinuit2/MinimumSeedGenerator.h"

namespace ROOT {

   namespace MyMinuit2 {


/** concrete implementation of the MinimumSeedGenerator interface; used within
    ModularFunctionMinimizer;
 */

class MnSeedGenerator : public MinimumSeedGenerator {

public:

  MnSeedGenerator() {}

  virtual ~MnSeedGenerator() {}

  virtual MinimumSeed operator()(const MnFcn&, const GradientCalculator&, const MnUserParameterState&, const MnStrategy&) const;

  virtual MinimumSeed operator()(const MnFcn&, const AnalyticalGradientCalculator&, const MnUserParameterState&, const MnStrategy&) const;

private:

};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_MnSeedGenerator
