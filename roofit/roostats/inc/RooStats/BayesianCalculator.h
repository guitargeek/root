// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_BayesianCalculator
#define ROOSTATS_BayesianCalculator

#include "TNamed.h"

#include "Math/IFunctionfwd.h"

#include "RooArgSet.h"

#include "RooStats/IntervalCalculator.h"

#include "RooStats/SimpleInterval.h"

class RooAbsData;
class RooAbsPdf;
class RooPlot;
class RooAbsReal;
class TF1;
class TH1;


namespace RooStats {

   class ModelConfig;
   class SimpleInterval;

   class BayesianCalculator : public IntervalCalculator, public TNamed {

   public:

      /// constructor
      BayesianCalculator( );

      /// constructor from a data set and a ModelConfig, which must outlive the calculator. The pdf, the prior pdf,
      /// the parameter of interest, the nuisance parameters, the constraint terms and the global observables are the
      /// ones of the ModelConfig.
      BayesianCalculator( RooAbsData& data,
                          ModelConfig& model );

      /// destructor
      ~BayesianCalculator() override;

      /// get the plot with option to get it normalized
      RooPlot* GetPosteriorPlot(bool norm = false, double precision = 0.01) const;

      /// return posterior pdf (object is managed by the user)
      RooAbsPdf* GetPosteriorPdf() const;
      /// return posterior function (object is managed by the BayesianCalculator class)
      RooAbsReal* GetPosteriorFunction() const;

      /// return the approximate posterior as histogram (TH1 object). Note the object is managed by the BayesianCalculator class
      TH1 * GetPosteriorHistogram() const;

      /// compute the interval. By Default a central interval is computed
      /// By using SetLeftTileFraction can control if central/ upper/lower interval
      /// For shortest interval use SetShortestInterval(true)
      SimpleInterval* GetInterval() const override ;

      void SetData( RooAbsData & data ) override {
         fData = &data;
         ClearAll();
      }


      /// set the model via the ModelConfig, which must outlive the calculator
      void SetModel( const ModelConfig& model ) override;

      /// set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      void SetTestSize( double size ) override {
         fSize = size;
         fValidInterval = false;
      }
      /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      void SetConfidenceLevel( double cl ) override { SetTestSize(1.-cl); }
      /// Get the size of the test (eg. rate of Type I error)
      double Size() const override { return fSize; }
      /// Get the Confidence level for the test
      double ConfidenceLevel() const override { return 1.-fSize; }

      /// set the fraction of probability content on the left tail
      /// Central limits use 0.5 (default case)
      /// for upper limits it is 0 and 1 for lower limit
      /// For shortest intervals a negative value (i.e. -1) must be given
      void SetLeftSideTailFraction(double leftSideFraction )  {fLeftSideFraction = leftSideFraction;}

      /// set the Bayesian calculator to compute the shortest interval (default is central interval)
      /// to switch off SetLeftSideTailFraction to the right value
      void SetShortestInterval() { fLeftSideFraction = -1; }

      /// set the precision of the Root Finder
      void SetBrfPrecision( double precision ) { fBrfPrecision = precision; }

      /// use directly the approximate posterior function obtained by binning it in nbins
      /// by default the cdf is used by integrating the posterior
      /// if a value of nbin <= 0 the cdf function will be used
      void SetScanOfPosterior(int nbin = 100) { fNScanBins = nbin; }

      /// set the number of iterations when running a MC integration algorithm
      /// If not set use default algorithmic values
      /// In case of ToyMC sampling of the nuisance the value is 100
      /// In case of using the GSL MCintegrations types the default value is
      /// defined in ROOT::Math::IntegratorMultiDimOptions::DefaultNCalls()
      virtual void SetNumIters(Int_t numIters)  { fNumIterations = numIters; }

      /// set the integration type (possible type are) :
      void SetIntegrationType(const char * type);

      /// return the mode (most probable value of the posterior function)
      double GetMode() const;

      // force the nuisance pdf when using the toy mc sampling
      void ForceNuisancePdf(RooAbsPdf & pdf) { fNuisancePdf = &pdf; }

   protected:

      void ClearAll() const;

      void ApproximatePosterior() const;

      void ComputeIntervalFromApproxPosterior(double lowerCutOff, double upperCutOff) const;

      void ComputeIntervalFromCdf(double lowerCutOff, double upperCutOff) const;

      void ComputeShortestInterval() const;

   private:

      // plan to replace the above: return a SimpleInterval integrating
      // over all other parameters except the one specified as argument
      // virtual SimpleInterval* GetInterval( RooRealVar* parameter  ) const { return 0; }

      // The content of the model.
      RooAbsPdf *Pdf() const;
      RooArgSet const &POI() const;
      RooAbsPdf *PriorPdf() const;
      RooArgSet NuisanceParameters() const;

      RooAbsData* fData;                         ///< data set
      ModelConfig const *fModel = nullptr;       ///<! the model (owned by its workspace)
      RooAbsPdf* fNuisancePdf;                   ///< nuisance pdf (needed when using nuisance sampling technique)

      mutable RooAbsPdf* fProductPdf;            ///< internal pointer to the product of the constraint terms
      mutable std::unique_ptr<RooAbsReal> fLogLike; ///< internal pointer to log likelihood function
      mutable RooAbsReal* fLikelihood;           ///< internal pointer to likelihood function
      mutable RooAbsReal* fIntegratedLikelihood; ///< integrated likelihood function, i.e - unnormalized posterior function
      mutable RooAbsPdf* fPosteriorPdf;          ///< normalized (on the poi) posterior pdf
      mutable ROOT::Math::IGenFunction * fPosteriorFunction;   ///< function representing the posterior
      mutable TF1 * fApproxPosterior;    ///< TF1 representing the scanned posterior function
      mutable double  fLower;          ///< computer lower interval bound
      mutable double  fUpper;          ///< upper interval bound
      mutable double  fNLLMin;         ///< minimum value of Nll
      double fSize;                      ///< size used for getting the interval
      double fLeftSideFraction;          ///< fraction of probability content on left side of interval
      double fBrfPrecision;              ///< root finder precision
      mutable int fNScanBins;            ///< number of bins to scan, if = -1 no scan is done (default)
      int fNumIterations;                ///< number of iterations (when using ToyMC)
      mutable bool    fValidInterval;

      TString fIntegrationType;

   protected:

      ClassDefOverride(BayesianCalculator,4)  // BayesianCalculator class

   };
}

#endif
