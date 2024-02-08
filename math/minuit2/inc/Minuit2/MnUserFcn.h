// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MnUserFcn
#define ROOT_MyMinuit2_MnUserFcn

#include "MyMinuit2/MnFcn.h"

namespace ROOT {

namespace MyMinuit2 {

class MnUserTransformation;

/**
   Wrapper used by Minuit of FCN interface
   containing a reference to the transformation object
 */
class MnUserFcn : public MnFcn {

public:
   MnUserFcn(const FCNBase &fcn, const MnUserTransformation &trafo, int ncall = 0)
      : MnFcn(fcn, ncall), fTransform(trafo)
   {
   }

   ~MnUserFcn() override {}

   double operator()(const MnAlgebraicVector &) const override;

private:
   const MnUserTransformation &fTransform;
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_MnUserFcn
