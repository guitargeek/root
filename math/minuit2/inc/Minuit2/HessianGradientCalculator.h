// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_HessianGradientCalculator
#define ROOT_MyMinuit2_HessianGradientCalculator

#include "MyMinuit2/GradientCalculator.h"
#include "MyMinuit2/MnMatrix.h"
#include <utility>

namespace ROOT {

namespace MyMinuit2 {

class MnFcn;
class MnUserTransformation;
class MnMachinePrecision;
class MnStrategy;

/**
   HessianGradientCalculator: class to calculate Gradient for Hessian
 */

class HessianGradientCalculator : public GradientCalculator {

public:
   HessianGradientCalculator(const MnFcn &fcn, const MnUserTransformation &par, const MnStrategy &stra)
      : fFcn(fcn), fTransformation(par), fStrategy(stra)
   {
   }

   ~HessianGradientCalculator() override {}

   FunctionGradient operator()(const MinimumParameters &) const override;

   FunctionGradient operator()(const MinimumParameters &, const FunctionGradient &) const override;

   std::pair<FunctionGradient, MnAlgebraicVector>
   DeltaGradient(const MinimumParameters &, const FunctionGradient &) const;

   const MnFcn &Fcn() const { return fFcn; }
   const MnUserTransformation &Trafo() const { return fTransformation; }
   const MnMachinePrecision &Precision() const;
   const MnStrategy &Strategy() const { return fStrategy; }

   unsigned int Ncycle() const;
   double StepTolerance() const;
   double GradTolerance() const;

private:
   const MnFcn &fFcn;
   const MnUserTransformation &fTransformation;
   const MnStrategy &fStrategy;
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_HessianGradientCalculator
