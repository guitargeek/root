/*
 * Project: RooFit
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef ROO_CONSTRAINT_SUM
#define ROO_CONSTRAINT_SUM

#include <RooAbsReal.h>
#include <RooListProxy.h>
#include <RooSetProxy.h>

class RooGaussian;
class RooPoisson;

class RooConstraintSum : public RooAbsReal {
public:
   RooConstraintSum() = default;
   RooConstraintSum(const char *name, const char *title, const RooArgSet &constraintSet, const RooArgSet &paramSet,
                    bool takeGlobalObservablesFromData = false);

   RooConstraintSum(const RooConstraintSum &other, const char *name = nullptr);
   TObject *clone(const char *newname) const override { return new RooConstraintSum(*this, newname); }

   const RooArgList &list() { return _set1; }

   bool setData(RooAbsData const &data, bool cloneData = true);
   /// \copydoc setData(RooAbsData const&, bool)
   bool setData(RooAbsData &data, bool cloneData = true) override
   {
      return setData(static_cast<RooAbsData const &>(data), cloneData);
   }

   void computeBatch(double *output, size_t size, RooFit::Detail::DataMap const &) const override;

   std::unique_ptr<RooAbsArg>
   compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext &ctx) const override;

   void translate(RooFit::Detail::CodeSquashContext &ctx) const override;

protected:
   double evaluate() const override;

private:

   bool hardcodeGaussian(RooGaussian const &gauss);
   bool hardcodePoisson(RooPoisson const &poiss);

   struct GaussianInfo {
      std::string obsName;
      double obsVal = 0.;
      double sigmaInvVal = 0.;
   };

   struct PoissonInfo {
      std::string obsName;
      double obsVal = 0.;
      bool noRounding = false;
   };

   std::vector<GaussianInfo> _gaussians;
   RooListProxy _gaussianParams;

   std::vector<PoissonInfo> _poissons;
   RooListProxy _poissonParams;

   RooListProxy _set1;                                ///< Set of constraint terms
   RooArgSet _paramSet;                               ///< Set of parameters to which constraints apply
   const bool _takeGlobalObservablesFromData = false; ///< If the global observable values are taken from data

   ClassDefOverride(RooConstraintSum, 0) // sum of -log of set of RooAbsPdf representing parameter constraints
};

#endif
