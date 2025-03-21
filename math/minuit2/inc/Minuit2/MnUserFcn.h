// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnUserFcn
#define ROOT_Minuit2_MnUserFcn

#include "Minuit2/MnFcn.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

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

   double operator()(const MnAlgebraicVector &) const override;

   // Access the parameter transformations.
   // For internal use in the Minuit2 implementation.
   const MnUserTransformation &transform() const { return fTransform; }

   double callWithTransformedParams(std::vector<double> const &vpar) const;

private:
   const MnUserTransformation &fTransform;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnUserFcn
