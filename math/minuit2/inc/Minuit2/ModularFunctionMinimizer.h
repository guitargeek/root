// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ModularFunctionMinimizer
#define ROOT_Minuit2_ModularFunctionMinimizer

#include "Minuit2/MnConfig.h"

#include <memory>
#include <vector>

namespace ROOT {

namespace Minuit2 {

class MinimumSeedGenerator;
class MinimumBuilder;
class MinimumSeed;
class MnFcn;
class GradientCalculator;
class MnUserParameterState;
class MnUserParameters;
class MnUserCovariance;
class MnStrategy;
class FCNBase;
class FCNGradientBase;
class FunctionMinimum;
class FumiliFCNBase;

//_____________________________________________________________
/**
   Base common class providing the API for all the minimizer
   Various Minimize methods are provided varying on the type of
   FCN function passesd and on the objects used for the parameters.
   The user may give FCN or FCN with Gradient,
   Parameter starting values and initial Error guess (sigma) (or "step size"),
   or Parameter starting values and initial covariance matrix;
   covariance matrix is stored in Upper triangular packed storage format,
   e.g. the Elements in the array are arranged like
   {a(0,0), a(0,1), a(1,1), a(0,2), a(1,2), a(2,2), ...},
   the size is nrow*(nrow+1)/2 (see also MnUserCovariance.h);
 */
class ModularFunctionMinimizer {

public:
   ModularFunctionMinimizer(std::unique_ptr<MinimumSeedGenerator> && minSeedGen, std::unique_ptr<MinimumBuilder> && minBuilder);

   virtual ~ModularFunctionMinimizer();

   FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &, const std::vector<double> &,
                            unsigned int stra = 1, unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &, const std::vector<double> &,
                            unsigned int stra = 1, unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &, unsigned int,
                            const std::vector<double> &, unsigned int stra = 1, unsigned int maxfcn = 0,
                            double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &, unsigned int,
                            const std::vector<double> &, unsigned int stra = 1, unsigned int maxfcn = 0,
                                    double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNBase &, const MnUserParameters &, const MnStrategy &,
                            unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameters &, const MnStrategy &,
                            unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNBase &, const MnUserParameters &, const MnUserCovariance &,
                            const MnStrategy &, unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameters &, const MnUserCovariance &,
                            const MnStrategy &, unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNBase &, const MnUserParameterState &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameterState &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   FunctionMinimum Minimize(const MnFcn &, const GradientCalculator &, const MinimumSeed &, const MnStrategy &,
                            unsigned int, double) const;

   const MinimumSeedGenerator &SeedGenerator() const { return *fMinSeedGen; }
   const MinimumBuilder &Builder() const { return *fMinBuilder; }
   MinimumBuilder &Builder() { return *fMinBuilder; }

private:
   std::unique_ptr<MinimumSeedGenerator> fMinSeedGen;
   std::unique_ptr<MinimumBuilder> fMinBuilder;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ModularFunctionMinimizer
