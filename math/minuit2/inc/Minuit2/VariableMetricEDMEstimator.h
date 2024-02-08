// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_VariableMetricEDMEstimator
#define ROOT_MyMinuit2_VariableMetricEDMEstimator

namespace ROOT {

namespace MyMinuit2 {

class FunctionGradient;
class MinimumError;

class VariableMetricEDMEstimator {

public:
   VariableMetricEDMEstimator() {}

   ~VariableMetricEDMEstimator() {}

   double Estimate(const FunctionGradient &, const MinimumError &) const;

private:
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_VariableMetricEDMEstimator
