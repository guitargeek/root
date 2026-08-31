/// \cond ROOFIT_INTERNAL

/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2023
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_FitHelpers_h
#define RooFit_FitHelpers_h

#include <memory>
#include <string>

class RooAbsData;
class RooDataHist;
class RooAbsPdf;
class RooAbsReal;
class RooCmdConfig;
class RooFitResult;
class RooLinkedList;

namespace RooFit {
namespace FitHelpers {

void defineMinimizationOptions(RooCmdConfig &pc);

std::unique_ptr<RooFitResult> minimize(RooAbsReal &model, RooAbsReal &nll, RooAbsData const &data, RooCmdConfig const &pc);

std::unique_ptr<RooAbsReal> createNLL(RooAbsPdf &pdf, RooAbsData &data, const RooLinkedList &cmdList);
std::unique_ptr<RooAbsReal> createChi2(RooAbsReal &real, RooDataHist &data, const RooLinkedList &cmdList);

std::unique_ptr<RooFitResult> fitTo(RooAbsReal &pdf, RooAbsData &data, const RooLinkedList &cmdList, bool chi2);

/// Build the evaluation-error message emitted when the chi-square denominator
/// of a bin is not positive, which makes the chi-square undefined there.
/// Shared by the legacy `RooChi2Var` and the vectorizing evaluation backends
/// so that the user gets the same actionable advice from both.
/// \param[in] binLabel Human-readable identification of the offending bin.
/// \param[in] nData Observed (weighted) event count in that bin.
/// \param[in] nPred Event count predicted by the model in that bin.
/// \param[in] expectedError Whether the denominator is the error predicted by
///            the model (Pearson) rather than the error of the data (Neyman).
///            The recommendation given to the user depends on it.
std::string chi2ZeroErrorBinMessage(std::string const &binLabel, double nData, double nPred, bool expectedError);

} // namespace FitHelpers
} // namespace RooFit

#endif

/// \endcond
