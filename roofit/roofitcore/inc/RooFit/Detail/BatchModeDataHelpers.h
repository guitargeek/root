/// \cond ROOFIT_INTERNAL

/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2022
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_BatchModeDataHelpers_h
#define RooFit_BatchModeDataHelpers_h

#include <RooFit/Detail/DataMap.h>

#include <ROOT/RSpan.hxx>
#include <string_view>

#include <functional>
#include <map>
#include <memory>
#include <stack>
#include <vector>

class RooAbsCategory;
class RooAbsData;
class RooSimultaneous;

namespace RooFit {
namespace Detail {

struct DataSpanInfo {
   std::span<const double> span;
   bool isGlobalObservable = false;
};

using DataSpanInfos = std::map<RooFit::Detail::DataKey, DataSpanInfo>;

DataSpanInfos getDataSpans(RooAbsData const &data, std::string const &rangeName, RooSimultaneous const *simPdf,
                           bool skipZeroWeights, std::stack<std::vector<double>> &buffers);

std::map<RooFit::Detail::DataKey, std::size_t>
determineOutputSizes(RooAbsArg const &topNode,
                     std::function<std::size_t(RooFit::Detail::DataKey)> const &inputSizeFunc);

} // namespace Detail
} // namespace RooFit

#endif

/// \endcond
