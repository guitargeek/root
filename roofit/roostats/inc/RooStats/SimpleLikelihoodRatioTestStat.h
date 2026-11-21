// @(#)root/roostats:$Id$
// Author: Kyle Cranmer and Sven Kreiss    June 2010
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_SimpleLikelihoodRatioTestStat
#define ROOSTATS_SimpleLikelihoodRatioTestStat

#include "Rtypes.h"

#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/TestStatistic.h"

namespace RooStats {

   class SimpleLikelihoodRatioTestStat : public TestStatistic {

   public:

      /// Constructor for proof, now removed. Do not use.
      SimpleLikelihoodRatioTestStat() = default;

      /// Takes null and alternate parameters from the current values of the
      /// variables of the pdfs. Can be overridden.
      ///
      /// The likelihoods are created from the pdfs of the models and their
      /// conditional observables. The constraint terms of the models are not
      /// included: they are the same for the null and the alternate
      /// hypothesis and cancel in the ratio. The models need to outlive the
      /// test statistic.
      SimpleLikelihoodRatioTestStat(const ModelConfig &nullModel, const ModelConfig &altModel)
         : fNullModel(&nullModel),
           fAltModel(&altModel)
      {
         std::unique_ptr<RooArgSet> allNullVars{fNullModel->GetPdf()->getVariables()};
         fNullParameters = allNullVars->snapshot();

         std::unique_ptr<RooArgSet> allAltVars{fAltModel->GetPdf()->getVariables()};
         fAltParameters = allAltVars->snapshot();
      }

      /// Takes null and alternate parameters from values in nullParameters
      /// and altParameters. Can be overridden.
      SimpleLikelihoodRatioTestStat(const ModelConfig &nullModel, const ModelConfig &altModel,
                                    const RooArgSet &nullParameters, const RooArgSet &altParameters)
         : fNullModel(&nullModel),
           fAltModel(&altModel),
           fNullParameters(nullParameters.snapshot()),
           fAltParameters(altParameters.snapshot())
      {
      }

      ~SimpleLikelihoodRatioTestStat() override {
         if (fNullParameters) delete fNullParameters;
         if (fAltParameters) delete fAltParameters;
      }

      static void SetAlwaysReuseNLL(bool flag);

      void SetReuseNLL(bool flag) { fReuseNll = flag ; }

      void SetNullParameters(const RooArgSet& nullParameters) {
         if (fNullParameters) delete fNullParameters;
         fFirstEval = true;
         fNullParameters = (RooArgSet*) nullParameters.snapshot();
      }

      void SetAltParameters(const RooArgSet& altParameters) {
         if (fAltParameters) delete fAltParameters;
         fFirstEval = true;
         fAltParameters = (RooArgSet*) altParameters.snapshot();
      }

      /// this should be possible with RooAbsCollection
      bool ParamsAreEqual() {
         if (!fNullParameters->equals(*fAltParameters)) return false;

         bool ret = true;

         for (auto nullIt = fNullParameters->begin(), altIt = fAltParameters->begin();
              nullIt != fNullParameters->end() && altIt != fAltParameters->end(); ++nullIt, ++altIt) {
            RooAbsReal *null = static_cast<RooAbsReal *>(*nullIt);
            RooAbsReal *alt = static_cast<RooAbsReal *>(*altIt);
            if (null->getVal() != alt->getVal())
               ret = false;
         }

         return ret;
      }


      double Evaluate(RooAbsData& data, RooArgSet& nullPOI) override;

      virtual void EnableDetailedOutput( bool e=true ) { fDetailedOutputEnabled = e; fDetailedOutput = nullptr; }
      const RooArgSet* GetDetailedOutput(void) const override { return fDetailedOutput.get(); }

      const TString GetVarName() const override {
         return "log(L(#mu_{1}) / L(#mu_{0}))";
      }

   private:

      static std::unique_ptr<RooAbsReal> createNLL(const ModelConfig &model, RooAbsData &data);

      const ModelConfig* fNullModel = nullptr;
      const ModelConfig* fAltModel = nullptr;
      RooArgSet* fNullParameters = nullptr;
      RooArgSet* fAltParameters = nullptr;
      bool fFirstEval = true;

      bool fDetailedOutputEnabled = false;
      std::unique_ptr<RooArgSet> fDetailedOutput; ///<!

      std::unique_ptr<RooAbsReal> fNllNull; ///<! transient copy of the null NLL
      std::unique_ptr<RooAbsReal> fNllAlt;  ///<!  transient copy of the alt NLL
      static bool fgAlwaysReuseNll ;
      bool fReuseNll = false;


   ClassDefOverride(SimpleLikelihoodRatioTestStat,5)
};

}

#endif
