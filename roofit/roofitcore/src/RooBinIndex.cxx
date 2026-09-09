/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2026
 *
 * Copyright (c) 2026, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file RooBinIndex.cxx
\class RooFit::Detail::RooBinIndex
\ingroup Roofitcore

A node in the computation graph that represents the flattened bin index
corresponding to the current values of a set of observables in a given binned
structure. See the doxygen comment of the class declaration for more
information.
**/

#include <RooFit/Detail/RooBinIndex.h>

#include <RooBinning.h>
#include <RooDataHist.h>
#include <RooFit/Detail/NormalizationHelpers.h>
#include <RooUniformBinning.h>

#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace RooFit::Detail {

RooBinIndex::RooBinIndex(const char *name, RooArgList const &vars, std::vector<RooAbsBinning const *> const &binnings,
                         std::vector<int> const &coefs)
   : RooAbsReal{name, name}, _vars{"vars", "vars", this}, _coefs{coefs}
{
   _vars.add(vars);
   _binnings.reserve(binnings.size());
   for (RooAbsBinning const *binning : binnings) {
      _binnings.emplace_back(binning->clone());
   }
}

RooBinIndex::RooBinIndex(const RooBinIndex &other, const char *name)
   : RooAbsReal{other, name}, _vars{"vars", this, other._vars}, _coefs{other._coefs}
{
   _binnings.reserve(other._binnings.size());
   for (auto const &binning : other._binnings) {
      _binnings.emplace_back(binning->clone());
   }
}

RooBinIndex *RooBinIndex::getOrCreate(CompileContext &ctx, RooAbsArg &owner, RooArgList const &vars,
                                      std::vector<RooAbsBinning const *> const &binnings,
                                      std::vector<int> const &coefs)
{
   // Build a key that encodes the content of the bin index calculation:
   // observable names, binning types and boundaries, and index coefficients.
   // Two consumers with equal keys can share the same node, even if they
   // reference different binning objects (e.g. the binning clones inside
   // different RooDataHists with identical structure).
   std::stringstream key;
   key.imbue(std::locale::classic());
   key << std::setprecision(std::numeric_limits<double>::max_digits10);

   std::string name = "binIndex" + std::to_string(ctx.sharedNodeCount());
   for (std::size_t i = 0; i < vars.size(); ++i) {
      RooAbsBinning const &binning = *binnings[i];
      key << vars[i].GetName() << ";" << binning.ClassName() << ";" << coefs[i] << ";";
      double const *boundaries = binning.array();
      for (int j = 0; j < binning.numBoundaries(); ++j) {
         key << boundaries[j] << ",";
      }
      key << ";";
      name += "_";
      name += vars[i].GetName();
   }

   if (RooAbsArg *existing = ctx.sharedNode(key.str())) {
      return static_cast<RooBinIndex *>(existing);
   }

   auto node = std::make_unique<RooBinIndex>(name.c_str(), vars, binnings, coefs);
   // The node is created fully compiled: its only servers are the observables
   // in the compiled computation graph.
   ctx.markAsCompiled(*node);
   ctx.registerSharedNode(key.str(), *node);
   RooBinIndex *out = node.get();
   owner.addOwnedComponents(std::move(node));
   return out;
}

RooBinIndex *RooBinIndex::getOrCreateForDataHist(CompileContext &ctx, RooAbsArg &owner, RooArgList const &vars,
                                                 RooDataHist const &dataHist)
{
   auto const &binnings = dataHist.getBinnings();

   // For now, shared bin index nodes are only created for one-dimensional
   // histograms. This covers the most important use case of HistFactory
   // models, where the observables of each channel are one-dimensional. For
   // more dimensions, the histogram-based classes don't agree on a common
   // flattened index ordering (the RooDataHist iterates the last variable
   // fastest, while e.g. the ParamHistFunc parameter list is ordered with the
   // first variable fastest), so sharing the index between different classes
   // would first require consolidating these conventions.
   if (vars.size() != 1 || binnings.size() != 1) {
      return nullptr;
   }

   std::vector<RooAbsBinning const *> binningPtrs;
   binningPtrs.reserve(binnings.size());
   for (std::size_t i = 0; i < binnings.size(); ++i) {
      RooAbsBinning const *binning = binnings[i].get();
      // Only static real-valued binnings are supported: a null binning stands
      // for a category dimension, and other binning implementations can
      // depend on external objects (like the parameters of a
      // RooParamBinning), which the bin index node would not track as
      // servers.
      if (!dynamic_cast<RooUniformBinning const *>(binning) && !dynamic_cast<RooBinning const *>(binning)) {
         return nullptr;
      }
      if (!dynamic_cast<RooAbsReal const *>(&vars[i])) {
         return nullptr;
      }
      binningPtrs.push_back(binning);
   }

   // Following the RooDataHist convention, the last variable is the
   // fastest-running index.
   std::vector<int> coefs(vars.size());
   int mult = 1;
   for (std::size_t i = vars.size(); i > 0; --i) {
      coefs[i - 1] = mult;
      mult *= binningPtrs[i - 1]->numBins();
   }

   return getOrCreate(ctx, owner, vars, binningPtrs, coefs);
}

double RooBinIndex::evaluate() const
{
   int idx = 0;
   for (std::size_t i = 0; i < _binnings.size(); ++i) {
      idx += _coefs[i] * _binnings[i]->binNumber(static_cast<RooAbsReal const &>(_vars[i]).getVal());
   }
   return idx;
}

void RooBinIndex::doEval(RooFit::EvalContext &ctx) const
{
   std::span<double> output = ctx.output();
   std::size_t n = output.size();

   _intBuffer.assign(n, 0);

   // Use the vectorized RooAbsBinning::binNumbers() to accumulate the
   // flattened bin index for each dimension, using the `coef` parameter to
   // multiply with the right index multiplication factor for each dimension.
   for (std::size_t i = 0; i < _binnings.size(); ++i) {
      _binnings[i]->binNumbers(ctx.at(&_vars[i]).data(), _intBuffer.data(), n, _coefs[i]);
   }

   for (std::size_t j = 0; j < n; ++j) {
      output[j] = _intBuffer[j];
   }
}

} // namespace RooFit::Detail
