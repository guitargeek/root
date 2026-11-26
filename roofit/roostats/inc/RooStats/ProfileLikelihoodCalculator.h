// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ProfileLikelihoodCalculator
#define ROOSTATS_ProfileLikelihoodCalculator

#include "RooStats/CombinedCalculator.h"

#include "RooStats/LikelihoodInterval.h"

namespace RooStats {

   class LikelihoodInterval;

 class ProfileLikelihoodCalculator : public CombinedCalculator {

   public:

      /// Default constructor (needed for I/O)
      ProfileLikelihoodCalculator();

      /// Constructor from data and a model configuration, which must outlive the calculator. The likelihood is the
      /// one of the ModelConfig, with its constraint terms. When using the calculator for performing an hypothesis
      /// test one needs to provide also a snapshot defining the null parameters and their value, see
      /// ModelConfig::SetSnapshot() and SetNullParameters(). There is no need to pass the alternate parameters. These
      /// will be obtained by the value maximizing the likelihood function.
      ProfileLikelihoodCalculator(RooAbsData& data, ModelConfig & model, double size = 0.05);


      ~ProfileLikelihoodCalculator() override;

      /// Return a likelihood interval. A global fit to the likelihood is performed and
      /// the interval is constructed using the profile likelihood ratio function of the POI.
      LikelihoodInterval* GetInterval() const override ;

      /// Return the hypothesis test result obtained from the likelihood ratio of the
      /// maximum likelihood value with the null parameters fixed to their values, with respect to keeping all parameters
      /// floating (global maximum likelihood value).
      HypoTestResult* GetHypoTest() const override;



   protected:

    /// clear internal fit result
    void DoReset() const;

    /// perform a global fit
    RooFit::OwningPtr<RooAbsReal> DoGlobalFit() const;

    /// minimize likelihood
    static RooFit::OwningPtr<RooFitResult> DoMinimizeNLL(RooAbsReal * nll);


    mutable std::unique_ptr<RooFitResult> fFitResult;  //// internal  result of global fit
    mutable bool fGlobalFitDone;        ///< flag to control if a global fit has been done


    ClassDefOverride(ProfileLikelihoodCalculator,3) // A concrete implementation of CombinedCalculator that uses the ProfileLikelihood ratio.

   };
}
#endif
