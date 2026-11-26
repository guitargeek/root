// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_CombinedCalculator
#define ROOSTATS_CombinedCalculator


#include "RooStats/IntervalCalculator.h"

#include "RooStats/HypoTestCalculator.h"

#include "RooStats/ModelConfig.h"

#include "RooAbsPdf.h"

#include "RooAbsData.h"

#include "RooArgSet.h"

// #ifndef ROO_WORKSPACE
// #include "RooWorkspace.h"
// #endif

namespace RooStats {

/** \class CombinedCalculator
   \ingroup Roostats

CombinedCalculator is an interface class for a tools which can produce both RooStats
HypoTestResults and ConfIntervals. The interface currently assumes that any such
calculator can be configured by specifying:

  - a model, given by a ModelConfig: the pdf, the parameters of interest, the nuisance
    parameters, the constraint terms and the global observables,
  - a data set,
  - a set of parameters of which specify the null (including values and const/non-const status),
  - a set of parameters of which specify the alternate (including values and const/non-const status).

The calculator does not copy the content of the ModelConfig: it keeps a pointer to it,
so the ModelConfig must outlive the calculator.

After configuring the calculator, one only needs to ask GetHypoTest() (which will
return a HypoTestResult pointer) or GetInterval() (which will return an ConfInterval pointer).

The concrete implementations of this interface should deal with the details of how
the nuisance parameters are dealt with (eg. integration vs. profiling) and which test-statistic is used (perhaps this should be added to the interface).

The motivation for this interface is that we hope to be able to specify the problem
in a common way for several concrete calculators.

*/


   class CombinedCalculator : public IntervalCalculator, public HypoTestCalculator {

   public:

      CombinedCalculator() = default;

      /// constructor from data and model configuration
      CombinedCalculator(RooAbsData &data, const ModelConfig &model, double size = 0.05) : fData(&data)
      {
         SetModel(model);
         SetTestSize(size);
      }

      /// Main interface to get a ConfInterval, pure virtual
      ConfInterval* GetInterval() const override = 0;

      /// main interface to get a HypoTestResult, pure virtual
      HypoTestResult* GetHypoTest() const override = 0;

      /// set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      void SetTestSize(double size) override {fSize = size;}

      /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      void SetConfidenceLevel(double cl) override {fSize = 1.-cl;}

      /// Get the size of the test (eg. rate of Type I error)
      double Size() const override {return fSize;}

      /// Get the Confidence level for the test
      double ConfidenceLevel()  const override {return 1.-fSize;}

      /// Set the DataSet, add to the workspace if not already there
      void SetData(RooAbsData & data) override {
         fData = &data;
      }

      /// Set the model, which must outlive the calculator. The parameter values of its snapshot are the null
      /// parameters.
      void SetModel(const ModelConfig & model) override {
         fModel = &model;
         if (model.GetSnapshot()) SetNullParameters(*model.GetSnapshot());
      }

      void SetNullModel( const ModelConfig &) override {  // to be understood what to do
      }
      void SetAlternateModel(const ModelConfig &) override {  // to be understood what to do
      }

      /// set parameter values for the null if using a common PDF
      virtual void SetNullParameters(const RooArgSet& set) {fNullParams.removeAll(); fNullParams.add(set);}

      /// set parameter values for the alternate if using a common PDF
      virtual void SetAlternateParameters(const RooArgSet& set) {fAlternateParams.removeAll(); fAlternateParams.add(set);}


   protected:

      ModelConfig const *GetModel() const { return fModel; }
      RooAbsPdf * GetPdf() const { return fModel ? fModel->GetPdf() : nullptr; }
      RooAbsData * GetData() const { return fData; }

      double fSize = 0.0; ///< size of the test (eg. specified rate of Type I error)

      ModelConfig const *fModel = nullptr; ///<! the model (owned by its workspace)
      RooAbsData * fData = nullptr;
      RooArgSet fNullParams;      ///< RooArgSet specifying null parameters for hypothesis test
      RooArgSet fAlternateParams; ///< RooArgSet specifying alternate parameters for hypothesis test


      ClassDefOverride(CombinedCalculator,3) // A base class that is for tools that can be both HypoTestCalculators and IntervalCalculators

   };
}


#endif
