/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_CONSTRAINT_SUM
#define ROO_CONSTRAINT_SUM

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooSetProxy.h"

class RooRealVar;
class RooArgList ;
class RooWorkspace ;

class RooGaussian;
class RooPoisson;

class RooConstraintSum : public RooAbsReal {
public:

  RooConstraintSum() {}
  RooConstraintSum(const char *name, const char *title, const RooArgSet& constraintSet, const RooArgSet& paramSet, bool takeGlobalObservablesFromData=false) ;

  RooConstraintSum(const RooConstraintSum& other, const char* name = nullptr);
  TObject* clone(const char* newname=nullptr) const override { return new RooConstraintSum(*this, newname); }

  const RooArgList& list() { return _set1 ; }

  bool setData(RooAbsData const& data, bool cloneData=true);
  /// \copydoc setData(RooAbsData const&, bool)
  bool setData(RooAbsData& data, bool cloneData=true) override {
    return setData(static_cast<RooAbsData const&>(data), cloneData);
  }

  void doEval(RooFit::EvalContext &) const override;

  std::unique_ptr<RooAbsArg> compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext & ctx) const override;

  struct GaussianInfo {
    double sigmaInvVal = 0.;
  };

  struct PoissonInfo {
    bool noRounding = false;
  };

  // Accessors for the hardcoded constraint terms, e.g. for code generation.
  std::vector<GaussianInfo> const &hardcodedGaussians() const { return _gaussians; }
  RooArgList const &gaussianParams() const { return _gaussianParams; }
  RooArgList const &gaussianObs() const { return _gaussianObs; }
  std::vector<PoissonInfo> const &hardcodedPoissons() const { return _poissons; }
  RooArgList const &poissonParams() const { return _poissonParams; }
  RooArgList const &poissonObs() const { return _poissonObs; }

protected:

  double evaluate() const override;

private:

  bool hardcodeGaussian(RooGaussian const &gauss);
  bool hardcodePoisson(RooPoisson const &poiss);

  // The observables are servers via the proxies and their values are read at
  // evaluation time: they are usually global observables, whose values can
  // change after the compilation for the fit (e.g. when they are taken from
  // the dataset, or in toy studies).
  std::vector<GaussianInfo> _gaussians;
  RooListProxy _gaussianParams;
  RooListProxy _gaussianObs;

  std::vector<PoissonInfo> _poissons;
  RooListProxy _poissonParams;
  RooListProxy _poissonObs;

  RooListProxy _set1 ;    ///< Set of constraint terms
  RooArgSet _paramSet ; ///< Set of parameters to which constraints apply
  const bool _takeGlobalObservablesFromData = false; ///< If the global observable values are taken from data

  ClassDefOverride(RooConstraintSum,0) // sum of -log of set of RooAbsPdf representing parameter constraints
};

#endif
