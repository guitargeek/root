/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file RooBinnedL.cxx
\class RooBinnedL
\ingroup Roofitcore

Implements a -log(likelihood) calculation from a dataset
(assumed to be binned) and a PDF. The NLL is calculated as
\f[
 \sum_\mathrm{data} -\log( \mathrm{pdf}(x_\mathrm{data}))
\f]
In extended mode, a
\f$ N_\mathrm{expect} - N_\mathrm{observed}*log(N_\mathrm{expect}) \f$ term is added.
**/

#include <RooFit/TestStatistics/RooBinnedL.h>

#include <RooFit/Detail/MathFuncs.h>

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsDataStore.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooChangeTracker.h"
#include "RooMsgService.h"

#include "TMath.h"

#include <stdexcept>
#include <string>

namespace RooFit {
namespace TestStatistics {

RooBinnedL::RooBinnedL(RooAbsPdf *pdf, RooAbsData *data, RooFit::ZeroPredictionMode zeroPredMode, double zeroPredDelta)
   : RooAbsL(RooAbsL::ClonePdfData{pdf, data}, data->numEntries(), 1),
     _zeroPredMode{zeroPredMode},
     _zeroPredDelta{zeroPredDelta}
{
   // pdf must be a RooRealSumPdf representing a yield vector for a binned likelihood calculation
   if (!dynamic_cast<RooRealSumPdf *>(pdf)) {
      throw std::logic_error("RooBinnedL can only be created from pdf of type RooRealSumPdf!");
   }

   // Retrieve and cache bin widths needed to convert unnormalized binned pdf values back to yields

   // The Active label will disable pdf integral calculations
   pdf->setAttribute("BinnedLikelihoodActive");

   RooArgSet params;
   pdf->getParameters(data->get(), params);
   paramTracker_ = std::make_unique<RooChangeTracker>("chtracker", "change tracker", params, true);

   std::unique_ptr<RooArgSet> obs(pdf->getObservables(data));
   if (obs->size() != 1) {
      throw std::logic_error(
         "RooBinnedL can only be created from combination of pdf and data which has exactly one observable!");
   } else {
      RooRealVar *var = static_cast<RooRealVar *>(obs->first());
      std::list<double> *boundaries = pdf->binBoundaries(*var, var->getMin(), var->getMax());
      std::list<double>::iterator biter = boundaries->begin();
      _binw.resize(boundaries->size() - 1);
      double lastBound = (*biter);
      ++biter;
      int ibin = 0;
      while (biter != boundaries->end()) {
         _binw[ibin] = (*biter) - lastBound;
         lastBound = (*biter);
         ibin++;
         ++biter;
      }
   }
}

RooBinnedL::~RooBinnedL() = default;

//////////////////////////////////////////////////////////////////////////////////
/// Calculate and return likelihood on subset of data from firstEvent to lastEvent
/// processed with a step size of 'stepSize'. If this an extended likelihood and
/// and the zero event is processed the extended term is added to the return
/// likelihood.
//
ROOT::Math::KahanSum<double>
RooBinnedL::evaluatePartition(Section bins, std::size_t /*components_begin*/, std::size_t /*components_end*/)
{
   // Throughout the calculation, we use Kahan's algorithm for summing to
   // prevent loss of precision - this is a factor four more expensive than
   // straight addition, but since evaluating the PDF is usually much more
   // expensive than that, we tolerate the additional cost...
   ROOT::Math::KahanSum<double> result;

   // Do not reevaluate likelihood if parameters nor event range have changed
   if (!paramTracker_->hasChanged(true) && bins == lastSection_ &&
       (cachedResult_.Sum() != 0 || cachedResult_.Carry() != 0))
      return cachedResult_;

   ROOT::Math::KahanSum<double> sumWeight;
   auto numEvalErrorsBefore = RooAbsReal::numEvalErrors();

   for (std::size_t i = bins.begin(N_events_); i < bins.end(N_events_); ++i) {

      data_->get(i);

      double eventWeight = data_->weight();

      // Calculate log(Poisson(N|mu) for this bin
      double N = eventWeight;
      double mu = pdf_->getVal() * _binw[i];

      if (_zeroPredMode != RooFit::ZeroPredictionMode::NaN) {
         // The zero-prediction regularization requested with
         // RooFit::ZeroPrediction() is active.
         if (N > 0. && mu < _zeroPredDelta) {
            if (_zeroPredMode == RooFit::ZeroPredictionMode::Error) {
               throw std::runtime_error(std::string{"RooBinnedL("} + pdf_->GetName() + "): observed " +
                                        std::to_string(N) + " events in bin " + std::to_string(i) +
                                        " where the model predicts " + std::to_string(mu) +
                                        " (ZeroPrediction(\"error\") is active)");
            }
            if (!_zeroPredWarned) {
               _zeroPredWarned = true;
               RooMsgService::instance().log(nullptr, RooFit::WARNING, RooFit::Eval)
                  << "RooBinnedL(" << pdf_->GetName() << "): bin " << i << " has a zero or tiny model prediction ("
                  << mu << ") with nonzero data (" << N
                  << "): the likelihood is regularized according to RooFit::ZeroPrediction (delta=" << _zeroPredDelta
                  << ")." << std::endl;
            }
         }

         double term = RooFit::Detail::MathFuncs::nllBinnedRegularized(mu, N, false, static_cast<int>(_zeroPredMode),
                                                                       _zeroPredDelta);
         sumWeight += eventWeight;
         result += term;
      } else if (mu <= 0 && N > 0) {

         // Catch error condition: data present where zero events are predicted
         RooAbsReal::logEvalError(nullptr, GetName().c_str(),
                                  TString::Format("Observed %f events in bin %zu with zero event yield", N, i));

      } else if (std::abs(mu) < 1e-10 && std::abs(N) < 1e-10) {

         // Special handling of this case since log(Poisson(0,0)=0 but can't be calculated with usual log-formula
         // since log(mu)=0. No update of result is required since term=0.

      } else {

         double term = -1 * (-mu + N * log(mu) - TMath::LnGamma(N + 1));

         sumWeight += eventWeight;
         result += term;
      }
   }

   // If part of simultaneous PDF normalize probability over
   // number of simultaneous PDFs: -sum(log(p/n)) = -sum(log(p)) + N*log(n)
   if (sim_count_ > 1) {
      result += sumWeight.Sum() * log(1.0 * sim_count_);
   }

   // At the end of the first full calculation, wire the caches
   if (_first) {
      _first = false;
      pdf_->wireAllCaches();
   }

   if ((RooAbsReal::evalErrorLoggingMode() == RooAbsReal::CollectErrors ||
        RooAbsReal::evalErrorLoggingMode() == RooAbsReal::CountErrors) &&
       numEvalErrorsBefore == RooAbsReal::numEvalErrors()) {
      cachedResult_ = result;
      lastSection_ = bins;
   }
   return result;
}

} // namespace TestStatistics
} // namespace RooFit
