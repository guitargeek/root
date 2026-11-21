// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additional Contributions: Giovanni Petrucciani
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ProfileLikelihoodTestStat
#define ROOSTATS_ProfileLikelihoodTestStat


#include "Rtypes.h"

#include "RooStats/TestStatistic.h"


#include "RooRealVar.h"

#include "Math/MinimizerOptions.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"


namespace RooStats {

  class ProfileLikelihoodTestStat : public TestStatistic{

     enum LimitType {twoSided, oneSided, oneSidedDiscovery};

   public:
      /// The likelihood is created with ModelConfig::createNLL(), such that
      /// it includes the constraint terms of the model and considers its
      /// global and conditional observables. The model needs to outlive the
      /// test statistic.
      ProfileLikelihoodTestStat(const ModelConfig &model)
         : fModel(&model),
           fStrategy(::ROOT::Math::MinimizerOptions::DefaultStrategy()),
           fTolerance(std::max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance())),
           fPrintLevel(::ROOT::Math::MinimizerOptions::DefaultPrintLevel()),
           fLOffset(RooStats::NLLOffsetMode())
      {
         // avoid default tolerance to be too small (1. is default in RooMinimizer)
      }

     ~ProfileLikelihoodTestStat() override {
       if(fCachedBestFitParams) delete fCachedBestFitParams;
     }

     void SetOneSided(bool flag=true) {fLimitType = (flag ? oneSided : twoSided);}
     void SetOneSidedDiscovery(bool flag=true) {fLimitType = (flag ? oneSidedDiscovery : twoSided);}
     void SetSigned(bool flag=true) {fSigned = flag;}  // +/- t_mu instead of t_mu>0 with one-sided settings

     bool IsTwoSided() const { return fLimitType == twoSided; }
     bool IsOneSidedDiscovery() const { return fLimitType == oneSidedDiscovery; }

     static void SetAlwaysReuseNLL(bool flag);

     void SetReuseNLL(bool flag) { fReuseNll = flag ; }
     void SetLOffset(bool flag=true) { fLOffset = flag ? "initial" : "none"; }
     void SetLOffset(std::string const &mode) { fLOffset = mode; }

     void SetMinimizer(const char* minimizer){ fMinimizer=minimizer;}
     void SetStrategy(Int_t strategy){fStrategy=strategy;}
     void SetTolerance(double tol){fTolerance=tol;}
     void SetPrintLevel(Int_t printlevel){fPrintLevel=printlevel;}

     /// Main interface to evaluate the test statistic on a dataset
     double Evaluate(RooAbsData& data, RooArgSet& paramsOfInterest) override {
        return EvaluateProfileLikelihood(0, data, paramsOfInterest);
     }

     /// evaluate  the profile likelihood ratio (type = 0) or the minimum of likelihood (type=1) or the conditional LL (type = 2)
     virtual double EvaluateProfileLikelihood(int type, RooAbsData &data, RooArgSet & paramsOfInterest);

     virtual void EnableDetailedOutput( bool e=true, bool withErrorsAndPulls=false ) {
        fDetailedOutputEnabled = e;
        fDetailedOutputWithErrorsAndPulls = withErrorsAndPulls;
        fDetailedOutput = nullptr;
     }
     /// Returns detailed output. The value returned by this function is updated after each call to Evaluate().
     /// The returned RooArgSet contains the following:
     ///
     ///  - the minimum nll, fitstatus and convergence quality for each fit
     ///  - for each fit and for each non-constant parameter, the value, error and pull of the parameter are stored
     ///
     const RooArgSet* GetDetailedOutput(void) const override {
      return fDetailedOutput.get();
     }

     virtual void SetVarName(const char* name) { fVarName = name; }
     const TString GetVarName() const override {return fVarName;}

     const ModelConfig &GetModel() const { return *fModel; }
     virtual RooAbsPdf * GetPdf() const { return fModel->GetPdf(); }

  private:

     std::unique_ptr<RooFitResult> GetMinNLL();
     bool minimizationNeeded(RooArgSet allParams) const;
     std::pair<double, int> minimizeNLL(std::string const &prefix);

      const ModelConfig* fModel = nullptr;
      std::unique_ptr<RooAbsReal> fNll; ///<!
      const RooArgSet* fCachedBestFitParams = nullptr;
      RooAbsData* fLastData = nullptr;
      LimitType fLimitType = twoSided;
      bool fSigned = false;

      /// this will store a snapshot of the unconditional nuisance
      /// parameter fit.
      bool fDetailedOutputEnabled = false;
      bool fDetailedOutputWithErrorsAndPulls = false;
      std::unique_ptr<RooArgSet> fDetailedOutput; ///<!

      TString fVarName = "Profile Likelihood Ratio";

      static bool fgAlwaysReuseNll ;
      bool fReuseNll = false;
      TString fMinimizer;
      Int_t fStrategy;
      double fTolerance;
      Int_t fPrintLevel;
      std::string fLOffset;

      ClassDefOverride(ProfileLikelihoodTestStat,0)   // implements the profile likelihood ratio as a test statistic to be used with several tools
   };
}


#endif
