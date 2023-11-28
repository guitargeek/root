/// \cond ROOFIT_INTERNAL

/*
 * Project: RooFit
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef ROO_BIN_INTEGRATOR
#define ROO_BIN_INTEGRATOR

#include "RooAbsIntegrator.h"
#include "RooNumIntConfig.h"
#include <vector>
#include <list>

class RooBinIntegrator : public RooAbsIntegrator {
public:
   RooBinIntegrator(const RooAbsFunc &function, int numBins = 100);
   RooBinIntegrator(const RooAbsFunc &function, const RooNumIntConfig &config);

   bool checkLimits() const override;
   double integral(const double *yvec = nullptr) override;

   using RooAbsIntegrator::setLimits;
   bool setLimits(double *xmin, double *xmax) override;
   bool setUseIntegrandLimits(bool flag) override
   {
      _useIntegrandLimits = flag;
      return true;
   }

private:
   friend class RooNumIntFactory;
   static void registerIntegrator(RooNumIntFactory &fact);
   RooBinIntegrator(const RooBinIntegrator &);

   // Numerical integrator workspace
   mutable std::vector<double> _xmin;      ///<! Lower integration bound
   mutable std::vector<double> _xmax;      ///<! Upper integration bound
   std::vector<std::vector<double>> _binb; ///<! list of bin boundaries
   int _numBins = 0;                       ///<! Size of integration range

   bool _useIntegrandLimits = false; ///< If true limits of function binding are ued

   std::vector<double> _x; ///<! do not persist
};

#endif

/// \endcond
