/// \cond ROOFIT_INTERNAL

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

#include "RooChi2Var.h"

#include "FitHelpers.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooCmdConfig.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooAbsCategory.h"
#include "RooAbsDataStore.h"

#include <ROOT/StringUtils.hxx>

#include <limits>
#include <ostream>
#include <sstream>
#include <string>

namespace {

/// Label a data-hist row by its index and, where available, by the values of
/// its observables, e.g. "bin 3 (x=0.175, cat=sig)".
std::string describeBin(std::size_t idx, RooArgSet const &row)
{
   std::stringstream ss;
   ss << "bin " << idx;
   std::stringstream coords;
   bool first = true;
   for (RooAbsArg *arg : row) {
      auto *real = dynamic_cast<RooAbsReal *>(arg);
      auto *cat = dynamic_cast<RooAbsCategory *>(arg);
      if (!real && !cat) {
         continue;
      }
      coords << (first ? "" : ", ") << arg->GetName() << "=";
      if (real) {
         coords << real->getVal();
      } else {
         coords << cat->getCurrentLabel();
      }
      first = false;
   }
   if (!first) {
      ss << " (" << coords.str() << ")";
   }
   return ss.str();
}

} // namespace

RooChi2Var::RooChi2Var(const char *name, const char *title, RooAbsReal &func, RooDataHist &data, bool extended,
                       RooDataHist::ErrorType etype, RooAbsTestStatistic::Configuration const &cfg)
   : RooAbsOptTestStatistic(name, title, func, data, RooArgSet{}, cfg),
     _etype{etype == RooAbsData::Auto ? (data.isNonPoissonWeighted() ? RooAbsData::SumW2 : RooAbsData::Expected)
                                      : etype},
     _funcMode{dynamic_cast<RooAbsPdf *>(&func) ? (extended ? ExtendedPdf : Pdf) : Function}
{
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooChi2Var::RooChi2Var(const RooChi2Var& other, const char* name) :
  RooAbsOptTestStatistic(other,name),
  _etype(other._etype),
  _funcMode(other._funcMode)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate chi^2 in partition from firstEvent to lastEvent using given stepSize
/// Throughout the calculation, we use Kahan's algorithm for summing to
/// prevent loss of precision - this is a factor four more expensive than
/// straight addition, but since evaluating the PDF is usually much more
/// expensive than that, we tolerate the additional cost...

double RooChi2Var::evaluatePartition(std::size_t firstEvent, std::size_t lastEvent, std::size_t stepSize) const
{
  double result(0);
  double carry(0);

  // With DataError(RooAbsData::None), all bin errors are zero by definition and
  // the chi-square is zero by convention (the vectorizing backends in
  // RooNLLVarNew::doEvalChi2() do the same).
  if (_etype == RooAbsData::None) {
    _evalCarry = 0.;
    return 0.;
  }

  // Set to true if at least one bin has a non-positive error, which makes the
  // chi-square undefined.
  bool hasUndefinedBin = false;

  // Also consider the composite case of multiple ranges
  std::vector<std::string> rangeTokens;
  if (!_rangeName.empty()) {
    rangeTokens = ROOT::Split(_rangeName, ",");
  }

  // Determine normalization factor depending on type of input function
  double normFactor(1) ;
  switch (_funcMode) {
  case Function: normFactor=1 ; break ;
  case Pdf: normFactor = _dataClone->sumEntries() ; break ;
  case ExtendedPdf: normFactor = (static_cast<RooAbsPdf*>(_funcClone))->expectedEvents(_dataClone->get()) ; break ;
  }

  // Loop over bins of dataset
  RooDataHist* hdata = static_cast<RooDataHist*>(_dataClone) ;
  for (auto i=firstEvent ; i<lastEvent ; i+=stepSize) {

    // get the data values for this event
    RooArgSet const *row = hdata->get(i);

    // Skip bins that are outside of the selected range
    bool doSelect(true) ;
    if (!_rangeName.empty()) {
      doSelect = false;
      // A row is selected if it is inside at least one complete named range.
      for (const auto &rangeName : rangeTokens) {
        bool inThisRange = true;
        for (const auto arg : *row) {
          if (!arg->inRange(rangeName.c_str())) {
            inThisRange = false;
            break;
          }
        }
        if (inThisRange) {
          doSelect = true;
          break;
        }
      }
    }
    if (!doSelect) continue ;

    const double nData = hdata->weight(i) ;

    const double nPdf = _funcClone->getVal(_normSet) * normFactor * hdata->binVolume(i) ;

    const double eExt = nPdf-nData ;


    double eInt ;
    if (_etype != RooAbsData::Expected) {
       double eIntLo;
       double eIntHi;
       hdata->weightError(eIntLo, eIntHi, _etype);
       eInt = (eExt > 0) ? eIntHi : eIntLo;
    } else {
      eInt = sqrt(nPdf) ;
    }

    // Skip cases where pdf=0 and there is no data
    if (0. == eInt * eInt && 0. == nData * nData && 0. == nPdf * nPdf) continue ;

    // The chi-square is undefined for a bin with zero error. Log an evaluation
    // error (so the minimizer's error handling kicks in) and make sure that the
    // final value is NaN instead of a plausible-looking number.
    if (0. == eInt * eInt) {
      const bool expectedError = (_etype == RooAbsData::Expected);
      const std::string msg =
        RooFit::FitHelpers::chi2ZeroErrorBinMessage(describeBin(i, *row), nData, nPdf, expectedError);
      logEvalError(msg.c_str());
      hasUndefinedBin = true;
      continue;
    }

//     std::cout << "Chi2Var[" << i << "] nData = " << nData << " nPdf = " << nPdf << " errorExt = " << eExt << " errorInt = " << eInt << " contrib = " << eExt*eExt/(eInt*eInt) << std::endl ;

    double term = eExt*eExt/(eInt*eInt) ;
    double y = term - carry;
    double t = result + y;
    carry = (t - result) - y;
    result = t;
  }

  _evalCarry = carry;
  if (hasUndefinedBin) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return result ;
}

/// \endcond
