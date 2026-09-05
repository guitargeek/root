/*
 * Project: RooFit
 * Authors:
 *   Garima Singh, CERN 2023
 *   Jonas Rembser, CERN 2023
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <RooFit/EvalContext.h>

#include <RooBatchCompute.h>
#include <RooRealVar.h>

#include <algorithm>
#include <stdexcept>

namespace {

// To avoid deleted move assignment.
template <class T>
void assignSpan(std::span<T> &to, std::span<T> const &from)
{
   to = from;
}

} // namespace

namespace RooFit {

std::span<const double> EvalContext::at(RooAbsArg const *arg, RooAbsArg const * /*caller*/)
{
   std::span<const double> out;

   if (!arg->hasDataToken()) {
      auto var = static_cast<RooRealVar const *>(arg);
      assignSpan(out, {&var->_value, 1});
   } else {
      std::size_t idx = arg->dataToken();
      out = _ctx[idx];
   }

   if (!_enableVectorBuffers || out.size() != 1) {
      return out;
   }

   if (_bufferIdx == _buffers.size()) {
      _buffers.emplace_back(RooBatchCompute::bufferSize);
   }

   double *buffer = _buffers[_bufferIdx].data();

   std::fill_n(buffer, RooBatchCompute::bufferSize, out[0]);
   assignSpan(out, {buffer, 1});

   ++_bufferIdx;

   return out;
}

void EvalContext::setConfig(RooAbsArg const *arg, RooBatchCompute::Config const &config)
{
   if (!arg->hasDataToken())
      return;
   std::size_t idx = arg->dataToken();
   _cfgs[idx] = config;
}

RooBatchCompute::Config EvalContext::config(RooAbsArg const *arg) const
{
   if (!arg->hasDataToken()) {
      return {};
   }
   std::size_t idx = arg->dataToken();
   return _cfgs[idx];
}

void EvalContext::resize(std::size_t n)
{
   _cfgs.resize(n);
   _ctx.resize(n);
   _supportRanges.resize(n, {0, std::numeric_limits<std::size_t>::max()});
}

/// \brief Declare that the values of the span registered for `arg` are
/// exactly zero outside of the index range [begin, end).
///
/// Consumers of the span can then skip the events outside of that range; see
/// e.g. RooAddPdf::doEval(). The declaration is reset whenever a new span is
/// registered for `arg`.
void EvalContext::setSupportRange(RooAbsArg const *arg, std::size_t begin, std::size_t end)
{
   if (!arg->hasDataToken())
      return;
   std::size_t idx = arg->dataToken();
   if (idx < _supportRanges.size()) {
      _supportRanges[idx] = {begin, end};
   }
}

/// \brief The index range outside of which the values of the span registered
/// for `arg` are known to be exactly zero, clamped to the span size.
///
/// For spans without a support declaration, this is simply the full range.
std::pair<std::size_t, std::size_t> EvalContext::supportRange(RooAbsArg const *arg) const
{
   if (!arg->hasDataToken() || arg->dataToken() >= _supportRanges.size()) {
      return {0, std::numeric_limits<std::size_t>::max()};
   }
   std::size_t idx = arg->dataToken();
   std::size_t size = _ctx[idx].size();
   return {std::min(_supportRanges[idx].first, size), std::min(_supportRanges[idx].second, size)};
}

/// \brief Sets the output value with an offset.
///
/// This function sets the output value with an offset for the given argument.
/// It should only be used in reducer nodes. Depending on the current
/// OffsetMode, the result will either be just the value, the value minus the
/// offset, of just the offset.
///
/// \param arg Pointer to the RooAbsArg object.
/// \param val The value to be set.
/// \param offset The offset value.
///
/// \throws std::runtime_error if the argument is not a reducer node.
void EvalContext::setOutputWithOffset(RooAbsArg const *arg, ROOT::Math::KahanSum<double> val,
                                      ROOT::Math::KahanSum<double> const &offset)
{
   if (!arg->isReducerNode()) {
      throw std::runtime_error("You can only use setOutputWithOffset() in reducer nodes!");
   }
   if (_offsetMode == OffsetMode::WithOffset) {
      val -= offset;
   } else if (_offsetMode == OffsetMode::OnlyOffset) {
      val = offset;
   }
   const_cast<double *>(_ctx[arg->dataToken()].data())[0] = val.Sum();
}

} // namespace RooFit
