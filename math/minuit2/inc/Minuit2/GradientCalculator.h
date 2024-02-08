// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_GradientCalculator
#define ROOT_MyMinuit2_GradientCalculator

#include "MyMinuit2/MnMatrixfwd.h"

namespace ROOT {

namespace MyMinuit2 {

class MinimumParameters;
class FunctionGradient;

/**
   interface class for gradient calculators
 */
class GradientCalculator {

public:
   virtual ~GradientCalculator() {}

   virtual FunctionGradient operator()(const MinimumParameters &) const = 0;

   virtual FunctionGradient operator()(const MinimumParameters &, const FunctionGradient &) const = 0;

   virtual bool Hessian(const MinimumParameters &, MnAlgebraicSymMatrix &) const { return false;}

   virtual bool G2(const MinimumParameters &, MnAlgebraicVector &) const { return false;}

};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_GradientCalculator
