// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_ScanMinimizer
#define ROOT_MyMinuit2_ScanMinimizer

#include "MyMinuit2/MnConfig.h"
#include "MyMinuit2/ModularFunctionMinimizer.h"
#include "MyMinuit2/ScanBuilder.h"
#include "MyMinuit2/SimplexSeedGenerator.h"

namespace ROOT {

namespace MyMinuit2 {

//_____________________________________________________________
/**
   Class implementing the required methods for a minimization using SCAN
   API is provided in the upper ROOT::MyMinuit2::ModularFunctionMinimizer class
 */

class ScanMinimizer : public ModularFunctionMinimizer {

public:
   ScanMinimizer() : fSeedGenerator(SimplexSeedGenerator()), fBuilder(ScanBuilder()) {}

   ~ScanMinimizer() override {}

   const MinimumSeedGenerator &SeedGenerator() const override { return fSeedGenerator; }
   const MinimumBuilder &Builder() const override { return fBuilder; }
   MinimumBuilder &Builder() override { return fBuilder; }

private:
   SimplexSeedGenerator fSeedGenerator;
   ScanBuilder fBuilder;
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_ScanMinimizer
