// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_InitialGradientCalculator
#define ROOT_MyMinuit2_InitialGradientCalculator

#include "MyMinuit2/GradientCalculator.h"

namespace ROOT {

   namespace MyMinuit2 {


class MnFcn;
class MnUserTransformation;
class MnMachinePrecision;
class MnStrategy;

/**
   Class to calculate an initial estimate of the gradient
 */
class InitialGradientCalculator : public GradientCalculator {

public:

  InitialGradientCalculator(const MnFcn& fcn, const MnUserTransformation& par,
                            const MnStrategy& stra) :
    fFcn(fcn), fTransformation(par), fStrategy(stra) {};

  virtual ~InitialGradientCalculator() {}

  virtual FunctionGradient operator()(const MinimumParameters&) const;

  virtual FunctionGradient operator()(const MinimumParameters&,
                                      const FunctionGradient&) const;

  const MnFcn& Fcn() const {return fFcn;}
  const MnUserTransformation& Trafo() const {return fTransformation;}
  const MnMachinePrecision& Precision() const;
  const MnStrategy& Strategy() const {return fStrategy;}

  unsigned int Ncycle() const;
  double StepTolerance() const;
  double GradTolerance() const;

private:

  const MnFcn& fFcn;
  const MnUserTransformation& fTransformation;
  const MnStrategy& fStrategy;
};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_InitialGradientCalculator
