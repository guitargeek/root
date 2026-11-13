/*
 * Project: RooFit
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_RooFitImplHelpers_h
#define RooFit_RooFitImplHelpers_h

#include <RooAbsArg.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooMsgService.h>

#include <RooNaNPacker.h>

#include <TMath.h>

#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

class RooAbsPdf;
class RooAbsData;

namespace RooHelpers {

std::pair<double, double> getRangeOrBinningInterval(RooAbsArg const *arg, const char *rangeName);

bool checkIfRangesOverlap(RooArgSet const &observables, std::vector<std::string> const &rangeNames);

std::string getColonSeparatedNameString(RooArgSet const &argSet, char delim = ':');
RooArgSet selectFromArgSet(RooArgSet const &, std::string const &names);

namespace Detail {

bool snapshotImpl(RooAbsCollection const &input, RooAbsCollection &output, bool deepCopy, RooArgSet const *observables);
RooAbsArg *cloneTreeWithSameParametersImpl(RooAbsArg const &arg, RooArgSet const *observables);

} // namespace Detail

/// Clone RooAbsArg object and reattach to original parameters.
template <class T>
std::unique_ptr<T> cloneTreeWithSameParameters(T const &arg, RooArgSet const *observables = nullptr)
{
   return std::unique_ptr<T>{static_cast<T *>(Detail::cloneTreeWithSameParametersImpl(arg, observables))};
}

std::string getRangeNameForSimComponent(std::string const &rangeName, bool splitRange, std::string const &catName);

struct BinnedLOutput {
   RooAbsPdf *binnedPdf = nullptr;
   bool isBinnedL = false;
};

BinnedLOutput getBinnedL(RooAbsPdf const &pdf);

void getSortedComputationGraph(RooAbsArg const &func, RooArgSet &out);

/// Check that `function` is constant (flat) inside each bin defined by the
/// sorted `boundaries` when scanning the observable `obs`. Several interior
/// points are sampled per bin and compared to the bin's first sample; if any
/// of them deviates by more than `relTol` (relative to the value scale), the
/// function is not flat and false is returned. The value of `obs` is restored
/// on return.
bool isFunctionFlatInBins(const RooAbsReal &function, RooAbsRealLValue &obs, std::span<const double> boundaries,
                          double relTol = 1e-9);

} // namespace RooHelpers

namespace RooFit::Detail {

std::string makeValidVarName(std::string const &in);

void replaceAll(std::string &inOut, std::string_view what, std::string_view with);

std::string makeSliceCutString(RooArgSet const &sliceDataSet);

// Inlined because this is called inside RooAbsPdf::getValV(), and therefore
// performance critical.
inline double normalizeWithNaNPacking(RooAbsPdf const &pdf, double rawVal, double normVal)
{

   if (normVal < 0. || (normVal == 0. && rawVal != 0)) {
      // Unreasonable normalisations. A zero integral can be tolerated if the function vanishes, though.
      const std::string msg = "p.d.f normalization integral is zero or negative: " + std::to_string(normVal);
      pdf.logEvalError(msg.c_str());
      return RooNaNPacker::packFloatIntoNaN(-normVal + (rawVal < 0. ? -rawVal : 0.));
   }

   if (rawVal < 0.) {
      std::stringstream ss;
      ss << "p.d.f value is less than zero (" << rawVal << "), trying to recover";
      pdf.logEvalError(ss.str().c_str());
      return RooNaNPacker::packFloatIntoNaN(-rawVal);
   }

   if (TMath::IsNaN(rawVal)) {
      pdf.logEvalError("p.d.f value is Not-a-Number");
      return rawVal;
   }

   return (rawVal == 0. && normVal == 0.) ? 0. : rawVal / normVal;
}

} // namespace RooFit::Detail

double toDouble(const char *s);
double toDouble(const std::string &s);

#endif
