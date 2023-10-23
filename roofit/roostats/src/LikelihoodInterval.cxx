// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*****************************************************************************
 * Project: RooStats
 * Package: RooFit/RooStats
 * @(#)root/roofit/roostats:$Id$
 * Authors:
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
 *
 *****************************************************************************/


/** \class RooStats::LikelihoodInterval
    \ingroup Roostats

   LikelihoodInterval is a concrete implementation of the RooStats::ConfInterval interface.
   It implements a connected N-dimensional intervals based on the contour of a likelihood ratio.
   The boundary of the interval is equivalent to a MINUIT/MINOS contour about the maximum likelihood estimator

   The interval does not need to be an ellipse (eg. it is not the HESSE error matrix).
   The level used to make the contour is the same as that used in MINOS, eg. it uses Wilks' theorem,
   which states that under certain regularity conditions the function -2* log (profile likelihood ratio) is asymptotically distributed as a chi^2 with N-dof, where
   N is the number of parameters of interest.


   Note, a boundary on the parameter space (eg. s>= 0) or a degeneracy (eg. mass of signal if Nsig = 0) can lead to violations of the conditions necessary for Wilks'
   theorem to be true.

   Also note, one can use any RooAbsReal as the function that will be used in the contour; however, the level of the contour
   is based on Wilks' theorem as stated above.


#### References

*  1. F. James., Minuit.Long writeup D506, CERN, 1998.

*/


#include "RooStats/LikelihoodInterval.h"
#include "RooStats/RooStatsUtils.h"

#include "RooAbsReal.h"
#include "RooMsgService.h"
#include "RooFitResult.h"

#include "Math/WrappedFunction.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/MinimizerOptions.h"
#include "RooFunctor.h"
#include "RooProfileLL.h"

#include "TMinuitMinimizer.h"

#include <string>
#include <algorithm>
#include <functional>
#include <cctype>   // need to use c version of toupper defined here


ClassImp(RooStats::LikelihoodInterval); ;

using namespace RooStats;
using namespace std;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor with name and title

LikelihoodInterval::LikelihoodInterval(const char* name) :
   ConfInterval(name), fBestFitParams(nullptr), fLikelihoodRatio(nullptr), fConfidenceLevel(0.95)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Alternate constructor taking a pointer to the profile likelihood ratio, parameter of interest and
/// optionally a snapshot of best parameter of interest for interval

LikelihoodInterval::LikelihoodInterval(const char* name, RooAbsReal* lr, const RooArgSet* params,  RooArgSet * bestParams) :
   ConfInterval(name),
   fParameters(*params),
   fBestFitParams(bestParams),
   fLikelihoodRatio(lr),
   fConfidenceLevel(0.95)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

LikelihoodInterval::~LikelihoodInterval()
{
   if (fBestFitParams) delete fBestFitParams;
   if (fLikelihoodRatio) delete fLikelihoodRatio;
}


////////////////////////////////////////////////////////////////////////////////
/// This is the main method to satisfy the RooStats::ConfInterval interface.
/// It returns true if the parameter point is in the interval.

bool LikelihoodInterval::IsInInterval(const RooArgSet &parameterPoint) const
{
   RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // Method to determine if a parameter point is in the interval
  if( !this->CheckParameters(parameterPoint) ) {
    std::cout << "parameters don't match" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }

  // make sure likelihood ratio is set
  if(!fLikelihoodRatio) {
    std::cout << "likelihood ratio not set" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }



  // set parameters
  SetParameters(&parameterPoint, std::unique_ptr<RooArgSet>{fLikelihoodRatio->getVariables()}.get());


  // evaluate likelihood ratio, see if it's bigger than threshold
  if (fLikelihoodRatio->getVal()<0){
    std::cout << "The likelihood ratio is < 0, indicates a bad minimum or numerical precision problems.  Will return true" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return true;
  }


  // here we use Wilks' theorem.
  if ( TMath::Prob( 2* fLikelihoodRatio->getVal(), parameterPoint.getSize()) < (1.-fConfidenceLevel) ){
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }


  RooMsgService::instance().setGlobalKillBelow(msglevel);

  return true;

}

////////////////////////////////////////////////////////////////////////////////
/// returns list of parameters

RooArgSet* LikelihoodInterval::GetParameters() const
{
   return new RooArgSet(fParameters);
}

////////////////////////////////////////////////////////////////////////////////
/// check that the parameters are correct

bool LikelihoodInterval::CheckParameters(const RooArgSet &parameterPoint) const
{
  if (parameterPoint.getSize() != fParameters.getSize() ) {
    std::cout << "size is wrong, parameters don't match" << std::endl;
    return false;
  }
  if ( ! parameterPoint.equals( fParameters ) ) {
    std::cout << "size is ok, but parameters don't match" << std::endl;
    return false;
  }
  return true;
}



////////////////////////////////////////////////////////////////////////////////
/// Compute lower limit, check first if limit has been computed
/// status is a boolean flag which will b set to false in case of error
/// and is true if calculation is successful
/// in case of error return also a lower limit value of zero

double LikelihoodInterval::LowerLimit(const RooRealVar& param, bool & status)
{
   double lower = 0;
   double upper = 0;
   status = FindLimits(param, lower, upper);
   return lower;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute upper limit, check first if limit has been computed
/// status is a boolean flag which will b set to false in case of error
/// and is true if calculation is successful
/// in case of error return also a lower limit value of zero

double LikelihoodInterval::UpperLimit(const RooRealVar& param, bool & status)
{
   double lower = 0;
   double upper = 0;
   status = FindLimits(param, lower, upper);
   return upper;
}


void LikelihoodInterval::ResetLimits() {
   // reset map with cached limits - called every time the test size or CL has been changed
   fLowerLimits.clear();
   fUpperLimits.clear();
}


bool LikelihoodInterval::CreateMinimizer() {
   // internal function to create minimizer object needed to find contours or interval limits
   // (running MINOS).
   // Minimizer must be Minuit or Minuit2

   RooProfileLL * profilell = dynamic_cast<RooProfileLL*>(fLikelihoodRatio);
   if (!profilell) return false;

   RooAbsReal & nll  = profilell->nll();
   // bind the nll function in the right interface for the Minimizer class
   // as a function of only the parameters (poi + nuisance parameters)

   std::unique_ptr<RooArgSet> partmp{profilell->getVariables()};
   // need to remove constant parameters
   RemoveConstantParameters(&*partmp);

   RooArgList params(*partmp);

   // need to restore values and errors for POI
   if (fBestFitParams) {
      for (int i = 0; i < params.getSize(); ++i) {
         RooRealVar & par =  (RooRealVar &) params[i];
         RooRealVar * fitPar =  (RooRealVar *) (fBestFitParams->find(par.GetName() ) );
         if (fitPar) {
            par.setVal( fitPar->getVal() );
            par.setError( fitPar->getError() );
         }
      }
   }


   fMinimizer = std::make_unique<RooMinimizer>(nll);
   fMinimizer->minimize(fMinimizer->minimizerType().c_str());
   if (std::unique_ptr<RooFitResult>(fMinimizer->save())->status() != 0) {
      ccoutE(Minimization) << "Error: Minimization failed  " << std::endl;
      return false;
   }
   return true;
}

bool LikelihoodInterval::FindLimits(const RooRealVar & param, double &lower, double & upper)
{
   // Method to find both lower and upper limits using MINOS
   // If cached values exist (limits have been already found) return them in that case
   // check first if limit has been computed
   // otherwise compute limit using MINOS
   // in case of failure lower and upper will maintain previous value (will not be modified)

   std::map<std::string, double>::const_iterator itrl = fLowerLimits.find(param.GetName());
   std::map<std::string, double>::const_iterator itru = fUpperLimits.find(param.GetName());
   if ( itrl != fLowerLimits.end() && itru != fUpperLimits.end() ) {
      lower = itrl->second;
      upper = itru->second;
      return true;
   }

   bool ret = true;
   if (!fMinimizer.get()) ret = CreateMinimizer();
   if (!ret) {
      ccoutE(Eval) << "Error returned from minimization of likelihood function - cannot find interval limits " << std::endl;
      return false;
   }

   // getting a 1D interval so ndf = 1
   double err_level = TMath::ChisquareQuantile(ConfidenceLevel(),1); // level for -2log LR
   err_level = err_level/2; // since we are using -log LR
   fMinimizer->setErrorLevel(err_level);

   ret = fMinimizer->minos(param);
   double elow = param.getAsymErrorLo();
   double eup = param.getAsymErrorHi();

   // WHEN error is zero normally is at limit
   if (elow == 0) {
      lower = param.getMin();
      ccoutW(Minimization) << "Warning: lower value for " << param.GetName() << " is at limit " << lower << std::endl;
   }
   else
      lower = param.getVal() + elow;  // elow is negative

   if (eup == 0) {
      ccoutW(Minimization) << "Warning: upper value for " << param.GetName() << " is at limit " << upper << std::endl;
      upper = param.getMax();
   }
   else
      upper = param.getVal() + eup;

   // store limits in the map
   // minos return error limit = minValue +/- error
   fLowerLimits[param.GetName()] = lower;
   fUpperLimits[param.GetName()] = upper;

   return true;
}


Int_t LikelihoodInterval::GetContourPoints(const RooRealVar & paramX, const RooRealVar & paramY, double * x, double *y, Int_t npoints ) {
   // use Minuit to find the contour of the likelihood function at the desired CL

   bool ret = true;
   if (!fMinimizer.get()) ret = CreateMinimizer();
   if (!ret) {
      coutE(Eval) << "LikelihoodInterval - Error returned creating minimizer for likelihood function - cannot find contour points " << std::endl;
      return 0;
   }

   std::unique_ptr<RooFitResult> res{fMinimizer->save()};
   const int ix = res->floatParsFinal().index(&paramX);
   const int iy = res->floatParsFinal().index(&paramY);

   ROOT::Math::Minimizer * minimizer = fMinimizer->fitter()->GetMinimizer();

   assert(minimizer);

   // getting a 2D contour so ndf = 2
   double cont_level = TMath::ChisquareQuantile(ConfidenceLevel(),2); // level for -2log LR
   cont_level = cont_level/2; // since we are using -log LR
   fMinimizer->setErrorLevel(cont_level);

   unsigned int ncp = npoints;
   coutI(Minimization)  << "LikelihoodInterval - Finding the contour of " << paramX.GetName() << " ( " << ix << " ) and " << paramY.GetName() << " ( " << iy << " ) " << std::endl;
   ret = minimizer->Contour(ix, iy, ncp, x, y );
   if (!ret) {
      coutE(Minimization) << "LikelihoodInterval - Error finding contour for parameters " << paramX.GetName() << " and " << paramY.GetName()  << std::endl;
      return 0;
   }
   if (int(ncp) < npoints) {
      coutW(Minimization) << "LikelihoodInterval -Warning - Less points calculated in contours np = " << ncp << " / " << npoints << std::endl;
   }

   return ncp;
 }
