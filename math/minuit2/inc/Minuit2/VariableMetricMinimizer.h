// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_VariableMetricMinimizer
#define ROOT_MyMinuit2_VariableMetricMinimizer

#include "MyMinuit2/MnConfig.h"
#include "MyMinuit2/ModularFunctionMinimizer.h"
#include "MyMinuit2/MnSeedGenerator.h"
#include "MyMinuit2/VariableMetricBuilder.h"

namespace ROOT {

   namespace MyMinuit2 {

      class BFGSMinimizerType {};

//______________________________________________________________________________
/**
    Instantiates the SeedGenerator and MinimumBuilder for
    Variable Metric Minimization method.
    API is provided in the upper ROOT::MyMinuit2::ModularFunctionMinimizer class

 */



class VariableMetricMinimizer : public ModularFunctionMinimizer {

   

public:

   class BFGSType {};

   VariableMetricMinimizer() : fMinSeedGen(MnSeedGenerator()),
                               fMinBuilder(VariableMetricBuilder()) {}

    VariableMetricMinimizer(BFGSType) :
       fMinSeedGen(MnSeedGenerator()),
       fMinBuilder(VariableMetricBuilder(VariableMetricBuilder::kBFGS)) {}

   ~VariableMetricMinimizer() {}

   const MinimumSeedGenerator& SeedGenerator() const {return fMinSeedGen;}
   const MinimumBuilder& Builder() const {return fMinBuilder;}
   MinimumBuilder& Builder()  {return fMinBuilder;}

private:

   MnSeedGenerator fMinSeedGen;
   VariableMetricBuilder fMinBuilder;
};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_VariableMetricMinimizer
