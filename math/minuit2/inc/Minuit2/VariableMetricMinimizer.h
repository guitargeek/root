// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_VariableMetricMinimizer
#define ROOT_Minuit2_VariableMetricMinimizer

#include "Minuit2/MnConfig.h"
#include "Minuit2/ModularFunctionMinimizer.h"
#include "Minuit2/MnSeedGenerator.h"
#include "Minuit2/VariableMetricBuilder.h"

namespace ROOT {

namespace Minuit2 {

class BFGSMinimizerType {
};

//______________________________________________________________________________
/**
    Instantiates the SeedGenerator and MinimumBuilder for
    Variable Metric Minimization method.
    API is provided in the upper ROOT::Minuit2::ModularFunctionMinimizer class

 */

class VariableMetricMinimizer : public ModularFunctionMinimizer {

public:
   class BFGSType {
   };

   VariableMetricMinimizer() : ModularFunctionMinimizer(std::make_unique<MnSeedGenerator>(), std::make_unique<VariableMetricBuilder>()) {}

   VariableMetricMinimizer(BFGSType)
      : ModularFunctionMinimizer(std::make_unique<MnSeedGenerator>(), std::make_unique<VariableMetricBuilder>(VariableMetricBuilder::kBFGS))
   {
   }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_VariableMetricMinimizer
