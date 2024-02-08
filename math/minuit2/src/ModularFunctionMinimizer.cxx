// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/ModularFunctionMinimizer.h"
#include "MyMinuit2/MinimumSeedGenerator.h"
#include "MyMinuit2/AnalyticalGradientCalculator.h"
#include "MyMinuit2/ExternalInternalGradientCalculator.h"
#include "MyMinuit2/Numerical2PGradientCalculator.h"
#include "MyMinuit2/MinimumBuilder.h"
#include "MyMinuit2/MinimumSeed.h"
#include "MyMinuit2/FunctionMinimum.h"
#include "MyMinuit2/MnUserParameterState.h"
#include "MyMinuit2/MnUserParameters.h"
#include "MyMinuit2/MnUserCovariance.h"
#include "MyMinuit2/MnUserTransformation.h"
#include "MyMinuit2/MnUserFcn.h"
#include "MyMinuit2/FCNBase.h"
#include "MyMinuit2/FCNGradientBase.h"
#include "MyMinuit2/MnStrategy.h"
#include "MyMinuit2/MnHesse.h"
#include "MyMinuit2/MnLineSearch.h"
#include "MyMinuit2/MnParabolaPoint.h"
#include "MyMinuit2/MnPrint.h"

namespace ROOT {

namespace MyMinuit2 {

// #include "MyMinuit2/MnUserParametersPrint.h"

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNBase &fcn, const std::vector<double> &par,
                                                   const std::vector<double> &err, unsigned int stra,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNBase and std::vector of double's for parameter values and errors (step sizes)
   MnUserParameterState st(par, err);
   MnStrategy strategy(stra);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNGradientBase &fcn, const std::vector<double> &par,
                                                   const std::vector<double> &err, unsigned int stra,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNGradientBase (use analytical gradient provided in FCN)
   // and std::vector of double's for parameter values and errors (step sizes)
   MnUserParameterState st(par, err);
   MnStrategy strategy(stra);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

// move nrow before cov to avoid ambiguities when using default parameters
FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNBase &fcn, const std::vector<double> &par,
                                                   unsigned int nrow, const std::vector<double> &cov, unsigned int stra,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNBase using std::vector for parameter error and
   // an std::vector of size n*(n+1)/2 for the covariance matrix  and n (rank of cov matrix)

   MnUserParameterState st(par, cov, nrow);
   MnStrategy strategy(stra);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNGradientBase &fcn, const std::vector<double> &par,
                                                   unsigned int nrow, const std::vector<double> &cov, unsigned int stra,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNGradientBase (use analytical gradient provided in FCN)
   // using std::vector for parameter error and
   // an std::vector of size n*(n+1)/2 for the covariance matrix  and n (rank of cov matrix)

   MnUserParameterState st(par, cov, nrow);
   MnStrategy strategy(stra);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNBase &fcn, const MnUserParameters &upar,
                                                   const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{
   // minimize from FCNBase and MnUserParameters object

   MnUserParameterState st(upar);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNGradientBase &fcn, const MnUserParameters &upar,
                                                   const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{
   // minimize from FCNGradientBase (use analytical gradient provided in FCN)  and MnUserParameters object

   MnUserParameterState st(upar);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNBase &fcn, const MnUserParameters &upar,
                                                   const MnUserCovariance &cov, const MnStrategy &strategy,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNBase and MnUserParameters and MnUserCovariance objects

   MnUserParameterState st(upar, cov);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNGradientBase &fcn, const MnUserParameters &upar,
                                                   const MnUserCovariance &cov, const MnStrategy &strategy,
                                                   unsigned int maxfcn, double toler) const
{
   // minimize from FCNGradientBase (use analytical gradient provided in FCN)  and
   // MnUserParameters MnUserCovariance objects

   MnUserParameterState st(upar, cov);
   return Minimize(fcn, st, strategy, maxfcn, toler);
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNBase &fcn, const MnUserParameterState &st,
                                                   const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{
   // minimize from a FCNBase and a MnUserparameterState - interface used by all the previous ones
   // based on FCNBase. Create in this case a NumericalGradient calculator
   // Create the minuit FCN wrapper (MnUserFcn) containing the transformation (int<->ext)

   // need MnUserFcn for difference int-ext parameters
   MnUserFcn mfcn(fcn, st.Trafo());
   Numerical2PGradientCalculator gc(mfcn, st.Trafo(), strategy);

   unsigned int npar = st.VariableParameters();
   if (maxfcn == 0)
      maxfcn = 200 + 100 * npar + 5 * npar * npar;
   MinimumSeed mnseeds = SeedGenerator()(mfcn, gc, st, strategy);

   return Minimize(mfcn, gc, mnseeds, strategy, maxfcn, toler);
}

// use Gradient here
FunctionMinimum ModularFunctionMinimizer::Minimize(const FCNGradientBase &fcn, const MnUserParameterState &st,
                                                   const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{
   // minimize from a FCNGradientBase and a MnUserParameterState -
   // interface based on FCNGradientBase (external/analytical gradients)
   // Create in this case an AnalyticalGradient calculator
   // Create the minuit FCN wrapper (MnUserFcn) containing the transformation (int<->ext)

   MnUserFcn mfcn(fcn, st.Trafo());
   std::unique_ptr<AnalyticalGradientCalculator> gc;
   if (fcn.gradParameterSpace() == GradientParameterSpace::Internal) {
        //        std::cout << "-- ModularFunctionMinimizer::Minimize: Internal parameter space" << std::endl;
        gc = std::unique_ptr<AnalyticalGradientCalculator>(new ExternalInternalGradientCalculator(fcn, st.Trafo()));
   } else {
        //        std::cout << "-- ModularFunctionMinimizer::Minimize: External parameter space" << std::endl;
        gc = std::make_unique<AnalyticalGradientCalculator>(fcn, st.Trafo());
   }

   unsigned int npar = st.VariableParameters();
   if (maxfcn == 0)
      maxfcn = 200 + 100 * npar + 5 * npar * npar;

   // compute seed (will use internally numerical gradient in case calculator does not implement g2 computations)
   MinimumSeed mnseeds = SeedGenerator()(mfcn, *gc, st, strategy);
   auto minimum = Minimize(mfcn, *gc, mnseeds, strategy, maxfcn, toler);

   return minimum;
}

FunctionMinimum ModularFunctionMinimizer::Minimize(const MnFcn &mfcn, const GradientCalculator &gc,
                                                   const MinimumSeed &seed, const MnStrategy &strategy,
                                                   unsigned int maxfcn, double toler) const
{
   // Interface used by all the others for the minimization using the base MinimumBuilder class
   // According to the contained type of MinimumBuilder the right type will be used

   MnPrint print("ModularFunctionMinimizer");

   const MinimumBuilder &mb = Builder();
   // std::cout << typeid(&mb).Name() << std::endl;
   double effective_toler = toler * mfcn.Up(); // scale tolerance with Up()
   // avoid tolerance too smalls (than limits)
   double eps = MnMachinePrecision().Eps2();
   if (effective_toler < eps)
      effective_toler = eps;

   // check if maxfcn is already exhausted
   // case already reached call limit
   if (mfcn.NumOfCalls() >= maxfcn) {
      print.Warn("Stop before iterating - call limit already exceeded");

      return FunctionMinimum(seed, std::vector<MinimumState>(1, seed.State()), mfcn.Up(),
                             FunctionMinimum::MnReachedCallLimit);
   }

   return mb.Minimum(mfcn, gc, seed, strategy, maxfcn, effective_toler);
}

} // namespace MyMinuit2

} // namespace ROOT
