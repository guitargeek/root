// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_AnalyticalGradientCalculator
#define ROOT_MyMinuit2_AnalyticalGradientCalculator

#include "MyMinuit2/GradientCalculator.h"

namespace ROOT {

   namespace MyMinuit2 {


class FCNGradientBase;
class MnUserTransformation;

class AnalyticalGradientCalculator : public GradientCalculator {

public:

  AnalyticalGradientCalculator(const FCNGradientBase& fcn, const MnUserTransformation& state) : fGradCalc(fcn), fTransformation(state) {}

  ~AnalyticalGradientCalculator() {}


  virtual FunctionGradient operator()(const MinimumParameters&) const;

  virtual FunctionGradient operator()(const MinimumParameters&,
                                      const FunctionGradient&) const;

  virtual bool CheckGradient() const;

private:

  const FCNGradientBase& fGradCalc;
  const MnUserTransformation& fTransformation;
};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_AnalyticalGradientCalculator
