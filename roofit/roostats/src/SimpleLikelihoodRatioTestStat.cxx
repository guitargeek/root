// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::SimpleLikelihoodRatioTestStat
    \ingroup Roostats

TestStatistic class that returns -log(L[null] / L[alt]) where
L is the likelihood.
It is often called as the LEP Test statistic.


*/

#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RooStatsUtils.h"

bool RooStats::SimpleLikelihoodRatioTestStat::fgAlwaysReuseNll = true ;

////////////////////////////////////////////////////////////////////////////////

void RooStats::SimpleLikelihoodRatioTestStat::SetAlwaysReuseNLL(bool flag) { fgAlwaysReuseNll = flag ; }

/// Create the likelihood of the pdf of the model, without the constraint
/// terms: they are the same for both hypotheses and cancel in the ratio.
std::unique_ptr<RooAbsReal>
RooStats::SimpleLikelihoodRatioTestStat::createNLL(const ModelConfig &model, RooAbsData &data)
{
   RooArgSet conditionalObs;
   if (model.GetConditionalObservables()) {
      conditionalObs.add(*model.GetConditionalObservables());
   }
   return std::unique_ptr<RooAbsReal>{
      model.GetPdf()->createNLL(data, RooFit::CloneData(false), RooFit::ConditionalObservables(conditionalObs))};
}

double RooStats::SimpleLikelihoodRatioTestStat::Evaluate(RooAbsData& data, RooArgSet& nullPOI) {

   if (fFirstEval && ParamsAreEqual()) {
      oocoutW(fNullParameters,InputArguments)
         << "Same RooArgSet used for null and alternate, so you must explicitly SetNullParameters and SetAlternateParameters or the likelihood ratio will always be 1."
         << std::endl;
   }

   fFirstEval = false;

   RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   bool reuse = (fReuseNll || fgAlwaysReuseNll) ;

   bool created = false ;
   if (!fNllNull) {
      fNllNull = createNLL(*fNullModel, data);
      created = true ;
   }
   if (reuse && !created) {
      fNllNull->setData(data, false) ;
   }

   // make sure we set the variables attached to this nll
   std::unique_ptr<RooArgSet> attachedSet{fNllNull->getVariables()};
   attachedSet->assign(*fNullParameters);
   attachedSet->assign(nullPOI);
   double nullNLL = fNllNull->getVal();

   //std::cout << std::endl << "SLRTS: null params:" << std::endl;
   //attachedSet->Print("v");


   if (!reuse) {
      fNllNull.reset();
   }

   created = false ;
   if (!fNllAlt) {
      fNllAlt = createNLL(*fAltModel, data);
      created = true ;
   }
   if (reuse && !created) {
      fNllAlt->setData(data, false) ;
   }
   // make sure we set the variables attached to this nll
   attachedSet = std::unique_ptr<RooArgSet>{fNllAlt->getVariables()};
   attachedSet->assign(*fAltParameters);
   double altNLL = fNllAlt->getVal();

   //std::cout << std::endl << "SLRTS: alt params:" << std::endl;
   //attachedSet->Print("v");


//   std::cout << std::endl << "SLRTS null NLL: " << nullNLL << "    alt NLL: " << altNLL << std::endl << std::endl;


   if (!reuse) {
      fNllAlt.reset();
   }



   // save this snapshot
   if( fDetailedOutputEnabled ) {
      if( !fDetailedOutput ) {
         fDetailedOutput = std::make_unique<RooArgSet>( *(new RooRealVar("nullNLL","null NLL",0)), "detailedOut_SLRTS" );
         fDetailedOutput->add( *(new RooRealVar("altNLL","alternate NLL",0)) );
      }
      fDetailedOutput->setRealValue( "nullNLL", nullNLL );
      fDetailedOutput->setRealValue( "altNLL", altNLL );
   }


   RooMsgService::instance().setGlobalKillBelow(msglevel);
   return nullNLL - altNLL;
}
