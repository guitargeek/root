// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_SimplexMinimizer
#define ROOT_MyMinuit2_SimplexMinimizer

#include "MyMinuit2/MnConfig.h"
#include "MyMinuit2/ModularFunctionMinimizer.h"
#include "MyMinuit2/SimplexBuilder.h"
#include "MyMinuit2/SimplexSeedGenerator.h"

#include <vector>

namespace ROOT {

   namespace MyMinuit2 {

//_____________________________________________________________
/**
   Class implementing the required methods for a minimization using Simplex.
   API is provided in the upper ROOT::MyMinuit2::ModularFunctionMinimizer class
 */

class SimplexMinimizer : public ModularFunctionMinimizer {

public:

   SimplexMinimizer() : fSeedGenerator(SimplexSeedGenerator()),
                        fBuilder(SimplexBuilder()) {}

   ~SimplexMinimizer() {}

   const MinimumSeedGenerator& SeedGenerator() const {return fSeedGenerator;}
   const MinimumBuilder& Builder() const {return fBuilder;}
   MinimumBuilder& Builder()  {return fBuilder;}

private:

   SimplexSeedGenerator fSeedGenerator;
   SimplexBuilder fBuilder;
};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_SimplexMinimizer
