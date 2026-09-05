/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2021
 *   Emmanouil Michalainas, CERN 2021
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file Evaluator.cxx
\class RooFit::Evaluator
\ingroup Roofitcore

Evaluates a RooAbsReal object in other ways than recursive graph
traversal. Currently, it is being used for evaluating a RooAbsReal object and
supplying the value to the minimizer, during a fit. The class scans the
dependencies and schedules the computations in a secure and efficient way. The
computations take place in the RooBatchCompute library and can be carried off
by either the CPU or a CUDA-supporting GPU. The Evaluator class takes care
of data transfers. An instance of this class is created every time
RooAbsPdf::fitTo() is called and gets destroyed when the fitting ends.
**/

#include <RooFit/Evaluator.h>

#include <RooAbsCategory.h>
#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooBatchCompute.h>
#include <RooMsgService.h>
#include <RooNameReg.h>
#include <RooSimultaneous.h>

#include <RooBatchCompute.h>

#include "BatchModeDataHelpers.h"
#include "RooFitImplHelpers.h"

#include <atomic>
#include <iomanip>
#include <numeric>
#include <unordered_set>

namespace RooFit {

namespace {

// To avoid deleted move assignment.
template <class T>
void assignSpan(std::span<T> &to, std::span<T> const &from)
{
   to = from;
}

void logArchitectureInfo(bool useGPU)
{
   // We have to exit early if the message stream is not active. Otherwise it's
   // possible that this function skips logging because it thinks it has
   // already logged, but actually it didn't.
   if (!RooMsgService::instance().isActive(nullptr, RooFit::Fitting, RooFit::INFO)) {
      return;
   }

   // Don't repeat logging architecture info if the useGPU option didn't change
   {
      // Second element of pair tracks whether this function has already been called
      static std::pair<bool, bool> lastUseGPU;
      if (lastUseGPU.second && lastUseGPU.first == useGPU)
         return;
      lastUseGPU = {useGPU, true};
   }

   auto log = [](std::string_view message) {
      oocxcoutI(static_cast<RooAbsArg *>(nullptr), Fitting) << message << std::endl;
   };

   if (RooBatchCompute::cpuArchitecture() == RooBatchCompute::Architecture::GENERIC) {
      log("using generic CPU library compiled with no vectorizations");
   } else {
      log(std::string("using CPU computation library compiled with -m") + RooBatchCompute::cpuArchitectureName());
   }
   if (useGPU) {
      log("using CUDA computation library");
   }
}

} // namespace

/// A struct used by the Evaluator to store information on the RooAbsArgs in
/// the computation graph.
struct NodeInfo {

   bool isScalar() const { return outputSize == 1; }

   RooAbsArg *absArg = nullptr;
   RooAbsArg::OperMode originalOperMode;

   std::shared_ptr<RooBatchCompute::AbsBuffer> buffer;
   std::size_t iNode = 0;
   int remClients = 0;
   bool copyAfterEvaluation = false;
   bool fromArrayInput = false;
   bool isVariable = false;
   bool isDirty = true;
   bool isCategory = false;
   bool hasLogged = false;
   bool computeInGPU = false;
   bool isValueServer = false; // if this node is a value server to the top node
   std::size_t outputSize = 1;

   // Range-restricted evaluation (see Evaluator::rangeRestrictionAnalysis()):
   // the node only computes the events [sliceBegin, sliceBegin + computeSize)
   // of the full event axis. For an unrestricted node, sliceBegin is zero and
   // computeSize equals outputSize.
   std::size_t sliceBegin = 0;
   std::size_t computeSize = 1;
   // Global event index that element 0 of the node's registered result span
   // corresponds to. Nonzero only for restricted nodes with compressed
   // buffers (a mask-gated product keeps a full-length buffer, so its frame
   // starts at zero).
   std::size_t frameBegin = 0;
   // A product gated by a data-only binary mask: it computes only its slice,
   // but its buffer stays full-length, with the mathematically exact value of
   // zero outside of the slice.
   bool isMaskedProduct = false;
   // The result span in the node's own frame, as last registered in the
   // evaluation context. Used to re-derive the per-client input spans when
   // restricted nodes are present.
   std::span<const double> canonicalSpan;
   // Change tracking (only maintained when restricted nodes are present): the
   // evaluation counter at which this node was last computed, and the global
   // event range in which its values changed in that computation.
   int lastComputeStamp = -1;
   std::size_t changedBegin = 0;
   std::size_t changedEnd = std::numeric_limits<std::size_t>::max();
   std::size_t lastSetValCount = std::numeric_limits<std::size_t>::max();
   int lastCatVal = std::numeric_limits<int>::max();
   double scalarBuffer = 0.0;
   std::vector<NodeInfo *> serverInfos;
   std::vector<NodeInfo *> clientInfos;

   /// Check the servers of a node that has been computed and release its
   /// resources if they are no longer needed. Buffers of nodes whose results
   /// are copied between host and device (copyAfterEvaluation) must not be
   /// released eagerly: their pinned host memory can still be the source of
   /// an asynchronous copy that was enqueued on the CUDA stream, and a new
   /// owner would overwrite it from the CPU without any stream ordering.
   /// Those buffers are released at the beginning of the next evaluation
   /// instead, after the stream was synchronized at the end of this one.
   void decrementRemainingClients()
   {
      if (--remClients == 0 && !fromArrayInput && !copyAfterEvaluation) {
         buffer.reset();
      }
   }
};

/// Construct a new Evaluator. The constructor analyzes and saves metadata about the graph,
/// useful for the evaluation of it that will be done later. In case the CUDA mode is selected,
/// there's also some CUDA-related initialization.
///
/// \param[in] absReal The RooAbsReal object that sits on top of the
///            computation graph that we want to evaluate.
/// \param[in] useGPU Whether the evaluation should be preferably done on the GPU.
Evaluator::Evaluator(const RooAbsReal &absReal, bool useGPU)
   : _topNode{const_cast<RooAbsReal &>(absReal)}, _useGPU{useGPU}
{
   RooBatchCompute::initCPU();
   if (useGPU && RooBatchCompute::initCUDA() != 0) {
      throw std::runtime_error("Can't create Evaluator in CUDA mode because RooBatchCompute CUDA could not be loaded!");
   }
   // Some checks and logging of used architectures
   logArchitectureInfo(_useGPU);

   _bufferManager = _useGPU ? RooBatchCompute::dispatchCUDA->createBufferManager()
                            : RooBatchCompute::dispatchCPU->createBufferManager();

   RooArgSet serverSet;
   ::RooHelpers::getSortedComputationGraph(_topNode, serverSet);

   _evalContextCPU.resize(serverSet.size());
   if (useGPU) {
      _evalContextCUDA.resize(serverSet.size());
   }

   std::map<RooFit::Detail::DataKey, NodeInfo *> nodeInfos;

   // Fill the ordered nodes list and initialize the node info structs.
   _nodes.reserve(serverSet.size());
   std::size_t iNode = 0;
   for (RooAbsArg *arg : serverSet) {

      _nodes.emplace_back();
      auto &nodeInfo = _nodes.back();
      _nodesMap[arg->namePtr()] = &nodeInfo;

      nodeInfo.absArg = arg;
      nodeInfo.originalOperMode = arg->operMode();
      nodeInfo.iNode = iNode;
      nodeInfos[arg] = &nodeInfo;

      if (dynamic_cast<RooRealVar const *>(arg)) {
         nodeInfo.isVariable = true;
      } else {
         arg->setDataToken(iNode);
      }
      if (dynamic_cast<RooAbsCategory const *>(arg)) {
         nodeInfo.isCategory = true;
      }

      ++iNode;
   }

   for (NodeInfo &info : _nodes) {
      info.serverInfos.reserve(info.absArg->servers().size());
      for (RooAbsArg *server : info.absArg->servers()) {
         if (server->isValueServer(*info.absArg)) {
            auto *serverInfo = nodeInfos.at(server);
            info.serverInfos.emplace_back(serverInfo);
            serverInfo->clientInfos.emplace_back(&info);
         }
      }
   }

   // Figure out which nodes are value servers to the top node
   _nodes.back().isValueServer = true; // the top node itself
   for (auto iter = _nodes.rbegin(); iter != _nodes.rend(); ++iter) {
      if (!iter->isValueServer)
         continue;
      for (auto &serverInfo : iter->serverInfos) {
         serverInfo->isValueServer = true;
      }
   }

   syncDataTokens();

   if (_useGPU) {
      // Create the single CUDA stream on which all GPU computations and data
      // transfers of this Evaluator are enqueued. The graph is evaluated in
      // topological order, so ordering the operations by the stream is enough
      // to guarantee correct results.
      _cudaStream = RooBatchCompute::dispatchCUDA->newCudaStream();
      RooBatchCompute::Config cfg;
      cfg.setCudaStream(_cudaStream);
      for (auto &info : _nodes) {
         _evalContextCUDA.setConfig(info.absArg, cfg);
      }
   }
}

/// If there are servers with the same name that got de-duplicated in the
/// `_nodes` list, we need to set their data tokens too. We find such nodes by
/// visiting the servers of every known node.
void Evaluator::syncDataTokens()
{
   for (NodeInfo &info : _nodes) {
      std::size_t iValueServer = 0;
      for (RooAbsArg *server : info.absArg->servers()) {
         if (server->isValueServer(*info.absArg)) {
            auto *knownServer = info.serverInfos[iValueServer]->absArg;
            if (knownServer->hasDataToken()) {
               server->setDataToken(knownServer->dataToken());
            }
            ++iValueServer;
         }
      }
   }
}

void Evaluator::setInput(std::string const &name, std::span<const double> inputArray, bool isOnDevice)
{
   if (isOnDevice && !_useGPU) {
      throw std::runtime_error("Evaluator can only take device array as input in CUDA mode!");
   }

   // Check if "name" is used in the computation graph. If yes, add the span to
   // the data map and set the node info accordingly.

   auto found = _nodesMap.find(RooNameReg::ptr(name.c_str()));

   if (found == _nodesMap.end())
      return;

   _needToUpdateOutputSizes = true;

   // Invalidate the caches that reducer nodes key on the input data, like the
   // cached sum of event weights in RooNLLVarNew. The counter is global so
   // that generation values can never alias between different Evaluators.
   {
      static std::atomic<std::size_t> nextInputGeneration{1};
      const std::size_t gen = ++nextInputGeneration;
      _evalContextCPU._inputGeneration = gen;
      _evalContextCUDA._inputGeneration = gen;
   }

   NodeInfo &info = *found->second;

   info.fromArrayInput = true;
   info.absArg->setDataToken(info.iNode);
   info.outputSize = inputArray.size();

   if (!_useGPU) {
      _evalContextCPU.set(info.absArg, inputArray);
      info.canonicalSpan = inputArray;
      info.frameBegin = 0;
      return;
   }

   if (info.outputSize <= 1) {
      // Empty or scalar observables from the data don't need to be
      // copied to the GPU.
      _evalContextCPU.set(info.absArg, inputArray);
      _evalContextCUDA.set(info.absArg, inputArray);
      return;
   }

   // For simplicity, we put the data on both host and device for
   // now. This could be optimized by inspecting the clients of the
   // variable.
   if (isOnDevice) {
      _evalContextCUDA.set(info.absArg, inputArray);
      auto gpuSpan = _evalContextCUDA.at(info.absArg);
      info.buffer = _bufferManager->makeCpuBuffer(gpuSpan.size());
      info.buffer->assignFromDevice(gpuSpan);
      _evalContextCPU.set(info.absArg, {info.buffer->hostReadPtr(), gpuSpan.size()});
   } else {
      _evalContextCPU.set(info.absArg, inputArray);
      auto cpuSpan = _evalContextCPU.at(info.absArg);
      info.buffer = _bufferManager->makeGpuBuffer(cpuSpan.size());
      info.buffer->assignFromHost(cpuSpan);
      _evalContextCUDA.set(info.absArg, {info.buffer->deviceReadPtr(), cpuSpan.size()});
   }
}

void Evaluator::updateOutputSizes()
{
   std::map<RooFit::Detail::DataKey, std::size_t> sizeMap;
   for (auto &info : _nodes) {
      if (info.fromArrayInput) {
         sizeMap[info.absArg] = info.outputSize;
      } else {
         // any buffer for temporary results is invalidated by resetting the output sizes
         info.buffer.reset();
      }
   }

   auto outputSizeMap =
      RooFit::BatchModeDataHelpers::determineOutputSizes(_topNode, [&](RooFit::Detail::DataKey key) -> int {
         auto found = sizeMap.find(key);
         return found != sizeMap.end() ? found->second : -1;
      });

   for (auto &info : _nodes) {
      info.outputSize = outputSizeMap.at(info.absArg);
      info.computeSize = info.outputSize;
      info.sliceBegin = 0;
      info.frameBegin = 0;
      info.isMaskedProduct = false;
      info.isDirty = true;
   }

   _hasRestrictedNodes = false;
   if (!_useGPU) {
      rangeRestrictionAnalysis();
   }
   _evalContextCPU._changeTracking = _hasRestrictedNodes;

   if (_useGPU) {
      markGPUNodes();
   }

   _needToUpdateOutputSizes = false;
}

/// Detect data-only binary-mask nodes (attribute "BinaryMask") whose set of
/// selected events forms one contiguous range of the event axis, and restrict
/// the evaluation of the products gated by them (attribute
/// "MaskGatedProduct"), and of the subgraphs that those products consume
/// exclusively, to that range. The masks depend only on the data, so the
/// analysis only has to run when new data is loaded and the restriction stays
/// valid for the whole fit. This keeps the total evaluation cost of e.g. a
/// simultaneous pdf compiled into a mixture (see
/// RooSimultaneous::compileForNormSet()) proportional to the number of
/// events, instead of the number of events times the number of channels.
void Evaluator::rangeRestrictionAnalysis()
{
   // Quick scan, so that computation graphs without masks don't pay anything.
   bool haveMasks = false;
   for (auto &info : _nodes) {
      haveMasks |= (info.outputSize > 1 && info.absArg->getAttribute("BinaryMask"));
   }
   if (!haveMasks)
      return;

   // Determine the contiguous event range that each data-only mask selects.
   std::map<NodeInfo const *, std::pair<std::size_t, std::size_t>> maskRanges;
   for (auto &info : _nodes) {
      if (info.outputSize <= 1 || info.fromArrayInput || !info.absArg->getAttribute("BinaryMask"))
         continue;
      bool dataOnly = true;
      for (NodeInfo *server : info.serverInfos) {
         dataOnly &= server->fromArrayInput;
      }
      if (!dataOnly)
         continue;

      // Evaluate the mask once, on the full range. It will not become dirty
      // again as long as the data doesn't change.
      computeCPUNode(info.absArg, info);
      info.isDirty = false;

      std::span<const double> vals = info.canonicalSpan;
      std::size_t begin = vals.size();
      std::size_t end = 0;
      for (std::size_t i = 0; i < vals.size(); ++i) {
         if (vals[i] != 0.0) {
            begin = std::min(begin, i);
            end = i + 1;
         }
      }
      // Only a non-empty contiguous range of selected events can be used. If
      // the selected events are scattered (e.g. data that is not sorted by
      // channel), the mask is not usable and the gated product is simply
      // evaluated on all events, like before.
      bool contiguous = begin < end;
      for (std::size_t i = begin; i < end && contiguous; ++i) {
         contiguous &= vals[i] != 0.0;
      }
      if (contiguous) {
         maskRanges[&info] = {begin, end};
      }
   }
   if (maskRanges.empty())
      return;

   // Fix the computed slice of each gated product to the range of its mask.
   // The product keeps a full-length buffer whose value outside of the slice
   // is exactly zero, because that's what multiplying with the mask yields.
   for (auto &info : _nodes) {
      if (info.outputSize <= 1 || !info.absArg->getAttribute("MaskGatedProduct"))
         continue;
      for (NodeInfo *server : info.serverInfos) {
         auto found = maskRanges.find(server);
         if (found != maskRanges.end()) {
            info.isMaskedProduct = true;
            info.sliceBegin = found->second.first;
            info.computeSize = found->second.second - found->second.first;
            break;
         }
      }
      _hasRestrictedNodes |= info.isMaskedProduct;
   }
   if (!_hasRestrictedNodes)
      return;

   // Propagate the demands down the graph, visiting clients before servers
   // (the node list is topologically sorted). Each vector node is then
   // restricted to the interval hull of what its clients read from it.
   std::vector<std::pair<std::size_t, std::size_t>> demand(_nodes.size(), {std::numeric_limits<std::size_t>::max(), 0});
   auto addDemand = [&](NodeInfo const *server, std::size_t begin, std::size_t end) {
      auto &d = demand[server->iNode];
      d.first = std::min(d.first, begin);
      d.second = std::max(d.second, end);
   };
   for (auto it = _nodes.rbegin(); it != _nodes.rend(); ++it) {
      NodeInfo &info = *it;
      // A scalar node (e.g. a reducer like the NLL class, or an integral)
      // consumes its vector inputs in full.
      if (info.outputSize == 1) {
         for (NodeInfo *server : info.serverInfos) {
            addDemand(server, 0, server->outputSize);
         }
         continue;
      }
      if (!info.isMaskedProduct) {
         auto const &d = demand[info.iNode];
         const bool hasDemand = d.first < d.second;
         const bool isTop = info.absArg == &_topNode;
         if (hasDemand && !isTop && !info.fromArrayInput) {
            const std::size_t begin = d.first;
            const std::size_t end = std::min(d.second, info.outputSize);
            if (end - begin < info.outputSize) {
               info.sliceBegin = begin;
               info.computeSize = end - begin;
               info.frameBegin = begin;
            }
         }
      }
      const std::size_t begin = info.sliceBegin;
      const std::size_t end = info.sliceBegin + info.computeSize;
      for (NodeInfo *server : info.serverInfos) {
         if (info.isMaskedProduct && maskRanges.find(server) != maskRanges.end()) {
            // The mask itself was already evaluated in full above and doesn't
            // depend on parameters, so it must not be restricted.
            addDemand(server, 0, server->outputSize);
         } else {
            addDemand(server, begin, end);
         }
      }
   }
}

/// When range-restricted nodes are present, the evaluation-context span of
/// every input has to be aligned to the frame of the node that is about to be
/// computed: element i of each input span must correspond to the global event
/// index `info.sliceBegin + i`. The canonical result spans of the servers
/// stay untouched, only the context entries are rewritten (they are rewritten
/// again before any other node is computed).
void Evaluator::prepareInputSpans(NodeInfo &info)
{
   // Scalar nodes (reducers, integrals) consume their vector inputs in full.
   const bool sliced = info.outputSize > 1;
   for (NodeInfo *server : info.serverInfos) {
      std::span<const double> const &canonical = server->canonicalSpan;
      if (canonical.size() <= 1 || !server->absArg->hasDataToken())
         continue;
      // The demand propagation in rangeRestrictionAnalysis() guarantees that
      // the frame of each server covers the slice of all of its clients. The
      // checks are cheap insurance against out-of-bounds reads in case that
      // invariant is ever broken.
      if (sliced && info.sliceBegin >= server->frameBegin &&
          info.sliceBegin - server->frameBegin + info.computeSize <= canonical.size()) {
         const std::size_t shift = info.sliceBegin - server->frameBegin;
         _evalContextCPU.set(server->absArg, {canonical.data() + shift, info.computeSize});
      } else {
         _evalContextCPU.set(server->absArg, canonical);
      }
      // Registering a span resets its support declaration, so it has to be
      // re-published. Only spans in the global frame carry it.
      if (server->isMaskedProduct && info.sliceBegin == server->frameBegin) {
         _evalContextCPU.setSupportRange(server->absArg, server->sliceBegin, server->sliceBegin + server->computeSize);
      }
   }
}

Evaluator::~Evaluator()
{
   for (auto &info : _nodes) {
      if (!info.isVariable) {
         info.absArg->resetDataToken();
      }
   }
   if (_cudaStream) {
      RooBatchCompute::dispatchCUDA->deleteCudaStream(_cudaStream);
   }
}

void Evaluator::computeCPUNode(const RooAbsArg *node, NodeInfo &info)
{
   using namespace Detail;

   const std::size_t nOut = info.outputSize;
   const std::size_t nCompute = info.computeSize; // equal to nOut unless the node is range-restricted

   double *buffer = nullptr;
   if (nOut == 1) {
      buffer = &info.scalarBuffer;
      if (_useGPU) {
         _evalContextCUDA.set(node, {buffer, nOut});
      }
   } else {
      if (!info.hasLogged && _useGPU) {
         RooAbsArg const &arg = *info.absArg;
         oocoutI(&arg, FastEvaluations) << "The argument " << arg.ClassName() << "::" << arg.GetName()
                                        << " could not be evaluated on the GPU because the class doesn't support it. "
                                           "Consider requesting or implementing it to benefit from a speed up."
                                        << std::endl;
         info.hasLogged = true;
      }
      if (!info.buffer) {
         const std::size_t nAlloc = info.isMaskedProduct ? nOut : nCompute;
         info.buffer = info.copyAfterEvaluation ? _bufferManager->makePinnedBuffer(nAlloc, _cudaStream)
                                                : _bufferManager->makeCpuBuffer(nAlloc);
         if (info.isMaskedProduct && nCompute != nOut) {
            // The value of a mask-gated product is exactly zero outside of
            // its slice: initialize the buffer once, only the slice is
            // rewritten by the evaluations.
            double *ptr = info.buffer->hostWritePtr();
            std::fill(ptr, ptr + nOut, 0.0);
         }
      }
      buffer = info.buffer->hostWritePtr();
   }
   // A mask-gated product keeps a full-length buffer and computes only its
   // slice; a compressed restricted node writes its slice at the beginning of
   // its (shorter) buffer.
   const std::size_t nRegister = info.isMaskedProduct ? nOut : nCompute;
   assignSpan(_evalContextCPU._currentOutput, {buffer + (info.isMaskedProduct ? info.sliceBegin : 0), nCompute});
   _evalContextCPU.set(node, {buffer, nRegister});
   if (info.isMaskedProduct) {
      _evalContextCPU.setSupportRange(node, info.sliceBegin, info.sliceBegin + info.computeSize);
   }
   assignSpan(info.canonicalSpan, {buffer, nRegister});
   if (nCompute > 1) {
      _evalContextCPU.enableVectorBuffers(true);
   }
   if (_hasRestrictedNodes) {
      prepareInputSpans(info);

      // Change tracking: determine in which global event range the values of
      // this node can change in this computation. For a node whose servers
      // were computed earlier in this evaluation, that's the union hull of
      // their changed ranges; otherwise (or for scalars) everything within
      // the node's own extent may change. A mask-gated product can only ever
      // change within its slice.
      if (nOut == 1) {
         info.changedBegin = 0;
         info.changedEnd = std::numeric_limits<std::size_t>::max();
      } else {
         std::size_t hullBegin = std::numeric_limits<std::size_t>::max();
         std::size_t hullEnd = 0;
         bool anyStampedServer = false;
         for (NodeInfo *server : info.serverInfos) {
            if (server->lastComputeStamp == _nEvaluations) {
               anyStampedServer = true;
               hullBegin = std::min(hullBegin, server->changedBegin);
               hullEnd = std::max(hullEnd, server->changedEnd);
            }
         }
         const std::size_t extentBegin = info.isMaskedProduct ? info.sliceBegin : info.frameBegin;
         const std::size_t extentEnd = extentBegin + nCompute;
         if (!anyStampedServer) {
            info.changedBegin = extentBegin;
            info.changedEnd = extentEnd;
         } else {
            info.changedBegin = std::max(extentBegin, hullBegin);
            info.changedEnd = std::min(extentEnd, hullEnd);
            if (info.changedBegin > info.changedEnd) {
               info.changedBegin = info.changedEnd = extentBegin;
            }
         }
         _evalContextCPU.setChangedRange(node, info.changedBegin, info.changedEnd);
      }
      info.lastComputeStamp = _nEvaluations;
   }
   if (info.isCategory) {
      auto nodeAbsCategory = static_cast<RooAbsCategory const *>(node);
      if (nOut == 1) {
         buffer[0] = nodeAbsCategory->getCurrentIndex();
      } else {
         throw std::runtime_error("RooFit::Evaluator - non-scalar category values are not supported!");
      }
   } else {
      auto nodeAbsReal = static_cast<RooAbsReal const *>(node);
      nodeAbsReal->doEval(_evalContextCPU);
   }
   _evalContextCPU.resetVectorBuffers();
   _evalContextCPU.enableVectorBuffers(false);
   if (info.copyAfterEvaluation) {
      // The deviceReadPtr() call triggers the copy of the result to the GPU.
      // The copy is ordered by the CUDA stream, so GPU clients enqueued later
      // will see the result without any further synchronization.
      _evalContextCUDA.set(node, {info.buffer->deviceReadPtr(), nOut});
   }
}

/// Process a variable in the computation graph. This is a separate non-inlined
/// function such that we can see in performance profiles how long this takes.
void Evaluator::processVariable(NodeInfo &nodeInfo)
{
   RooAbsArg *node = nodeInfo.absArg;
   auto *var = static_cast<RooRealVar const *>(node);
   if (nodeInfo.lastSetValCount != var->valueResetCounter()) {
      nodeInfo.lastSetValCount = var->valueResetCounter();
      for (NodeInfo *clientInfo : nodeInfo.clientInfos) {
         clientInfo->isDirty = true;
      }
      computeCPUNode(node, nodeInfo);
      nodeInfo.isDirty = false;
   }
}

/// Process a category in the computation graph. This is a separate non-inlined
/// function such that we can see in performance profiles how long this takes.
void Evaluator::processCategory(NodeInfo &nodeInfo)
{
   RooAbsArg *node = nodeInfo.absArg;
   auto *cat = static_cast<RooAbsCategory const *>(node);
   if (nodeInfo.lastCatVal != cat->getCurrentIndex()) {
      nodeInfo.lastCatVal = cat->getCurrentIndex();
      for (NodeInfo *clientInfo : nodeInfo.clientInfos) {
         clientInfo->isDirty = true;
      }
      computeCPUNode(node, nodeInfo);
      nodeInfo.isDirty = false;
   }
}

/// Flags all the clients of a given node dirty. This is a separate non-inlined
/// function such that we can see in performance profiles how long this takes.
void Evaluator::setClientsDirty(NodeInfo &nodeInfo)
{
   for (NodeInfo *clientInfo : nodeInfo.clientInfos) {
      clientInfo->isDirty = true;
   }
}

/// Returns the value of the top node in the computation graph
std::span<const double> Evaluator::run()
{
   if (_needToUpdateOutputSizes)
      updateOutputSizes();

   ++_nEvaluations;

   // Discard leftover deferred actions in case a previous evaluation was
   // aborted by an exception.
   _evalContextCPU._deferredActions.clear();
   _evalContextCUDA._deferredActions.clear();

   if (_useGPU) {
      return getValHeterogeneous();
   }

   for (auto &nodeInfo : _nodes) {
      if (!nodeInfo.fromArrayInput) {
         if (nodeInfo.isVariable) {
            processVariable(nodeInfo);
         } else if (nodeInfo.isCategory) {
            processCategory(nodeInfo);
         } else {
            if (nodeInfo.isDirty) {
               setClientsDirty(nodeInfo);
               computeCPUNode(nodeInfo.absArg, nodeInfo);
               nodeInfo.isDirty = false;
            }
         }
      }
   }

   for (auto &action : _evalContextCPU._deferredActions) {
      action();
   }
   _evalContextCPU._deferredActions.clear();

   // return the final output
   return _evalContextCPU.at(&_topNode);
}

/// Returns the value of the top node in the computation graph
std::span<const double> Evaluator::getValHeterogeneous()
{
   for (auto &info : _nodes) {
      info.remClients = info.clientInfos.size();
      if (info.buffer && !info.fromArrayInput) {
         info.buffer.reset();
      }
   }

   // Iterate over the nodes in topological order. Nodes that are computed on
   // the GPU only enqueue their computation on the single CUDA stream and
   // return immediately, so independent CPU nodes that come later in the
   // ordering naturally overlap with the GPU computations. Ordering by the
   // stream guarantees that GPU nodes see the results of their GPU servers,
   // and host-side reads of GPU results synchronize on the stream in the
   // buffer implementation.
   try {
      for (auto &info : _nodes) {
         if (!info.fromArrayInput) {
            if (info.computeInGPU) {
               assignToGPU(info);
            } else {
               computeCPUNode(info.absArg, info);
            }
         }

         // Release the buffers of server nodes that are no longer needed. For
         // device-only buffers this is safe to do right away even if GPU work
         // is still in flight, because any reuse of a released device buffer
         // happens through operations that are enqueued later on the same
         // stream. Pinned buffers are exempted from the eager release, see
         // the comment in NodeInfo::decrementRemainingClients().
         for (auto *serverInfo : info.serverInfos) {
            serverInfo->decrementRemainingClients();
         }
      }
   } catch (...) {
      // The evaluation was aborted, but readbacks that compute() calls
      // deferred may still be armed. Deliver them now, while the destination
      // memory in the nodes of the computation graph is guaranteed to be
      // alive, so that no armed readback survives into a later evaluation.
      try {
         RooBatchCompute::dispatchCUDA->synchronizeCudaStream(_cudaStream);
      } catch (...) {
         // The stream is in an unrecoverable error state. The deferred
         // readbacks are dropped together with the scratch memory when the
         // stream gets deleted.
      }
      _evalContextCUDA._deferredActions.clear();
      _evalContextCPU._deferredActions.clear();
      throw;
   }

   // Ensure that all enqueued GPU work has completed when run() returns. For
   // the usual likelihood evaluations this is mostly a no-op, because the
   // final reduction has synchronized the stream already. It also guarantees
   // that recycling the buffers at the beginning of the next evaluation is
   // safe, and it delivers the deferred readbacks like the evaluation error
   // counters.
   RooBatchCompute::dispatchCUDA->synchronizeCudaStream(_cudaStream);

   // Run the deferred actions now that all results have arrived on the host,
   // e.g. the logging of evaluation errors that were counted on the GPU.
   // Nodes evaluated on the CPU register their actions in the CPU context,
   // so both contexts are drained.
   for (auto *ctx : {&_evalContextCUDA, &_evalContextCPU}) {
      for (auto &action : ctx->_deferredActions) {
         action();
      }
      ctx->_deferredActions.clear();
   }

   // return the final value
   return _evalContextCUDA.at(&_topNode);
}

/// Enqueue the computation of a node on the GPU.
void Evaluator::assignToGPU(NodeInfo &info)
{
   using namespace Detail;

   auto node = static_cast<RooAbsReal const *>(info.absArg);

   const std::size_t nOut = info.outputSize;

   double *buffer = nullptr;
   if (nOut == 1) {
      buffer = &info.scalarBuffer;
      _evalContextCPU.set(node, {buffer, nOut});
   } else {
      info.buffer = info.copyAfterEvaluation ? _bufferManager->makePinnedBuffer(nOut, _cudaStream)
                                             : _bufferManager->makeGpuBuffer(nOut);
      buffer = info.buffer->deviceWritePtr();
   }
   assignSpan(_evalContextCUDA._currentOutput, {buffer, nOut});
   _evalContextCUDA.set(node, {buffer, nOut});
   node->doEval(_evalContextCUDA);
   if (info.copyAfterEvaluation) {
      // The hostReadPtr() call triggers the copy of the result to the host,
      // which waits for the enqueued computation via the CUDA stream.
      _evalContextCPU.set(node, {info.buffer->hostReadPtr(), nOut});
   }
}

/// Decides which nodes are assigned to the GPU in a CUDA fit.
void Evaluator::markGPUNodes()
{
   // Decide which nodes get evaluated on the GPU: we select nodes that support
   // CUDA evaluation and have at least one input of size greater than one.
   for (auto &info : _nodes) {
      info.computeInGPU = false;
      if (!info.absArg->canComputeBatchWithCuda()) {
         continue;
      }
      for (NodeInfo const *serverInfo : info.serverInfos) {
         if (serverInfo->outputSize > 1) {
            info.computeInGPU = true;
            break;
         }
      }
   }

   // In a second pass, figure out which nodes need to copy over their results.
   for (auto &info : _nodes) {
      info.copyAfterEvaluation = false;
      // scalar nodes don't need copying
      if (!info.isScalar()) {
         for (auto *clientInfo : info.clientInfos) {
            if (info.computeInGPU != clientInfo->computeInGPU) {
               info.copyAfterEvaluation = true;
               break;
            }
         }
      }
   }
}

/// \brief Sets the number of threads to use for the evaluation of a single node.
///
/// With a value greater than one, the computation functions and reductions of
/// the CPU backend process large batches multi-threaded, using up to the
/// given number of threads. Nodes evaluated on the CPU with fewer events than
/// an internal threshold are still evaluated single-threaded, so requesting
/// multiple threads never introduces scheduling overhead for small fits.
void Evaluator::setNThreads(int nThreads)
{
   for (auto &info : _nodes) {
      if (info.isVariable) {
         continue;
      }
      RooBatchCompute::Config cfg = _evalContextCPU.config(info.absArg);
      cfg.setNThreads(nThreads);
      _evalContextCPU.setConfig(info.absArg, cfg);
   }
}

/// Temporarily change the operation mode of a RooAbsArg until the
/// Evaluator gets deleted.
void Evaluator::setOperMode(RooAbsArg *arg, RooAbsArg::OperMode opMode)
{
   if (!_operModeChanges)
      _operModeChanges = std::make_unique<ChangeOperModeRAII>();
   _operModeChanges->change(arg, opMode);
}

// Change the operation modes of all RooAbsArgs in the computation graph.
// The changes are reset when the returned RAII object goes out of scope.
//
// We also walk transitively through value clients of the nodes to cover any
// node that RooAbsReal::doEval (the fallback scalar implementation) might
// inadvertently propagate the ADirty mode to via its recursive restore: that
// helper sets servers temporarily to AClean and then calls
// setOperMode(oldOperMode) to restore, which recurses to value clients when
// oldOperMode is ADirty. If we did not protect those clients here, any node
// outside the computation graph that shares a fundamental (e.g. a parameter
// like a RooRealVar) would be left permanently in ADirty after the first
// minimization, dramatically slowing down later scalar evaluations (for
// example on pdfs held by the legacy test statistics' internal cache).
std::unique_ptr<ChangeOperModeRAII> Evaluator::setOperModes(RooAbsArg::OperMode opMode)
{
   auto out = std::make_unique<ChangeOperModeRAII>();
   std::unordered_set<RooAbsArg *> visited;

   std::vector<RooAbsArg *> queue;
   queue.reserve(_nodes.size());
   for (auto &info : _nodes) {
      queue.push_back(info.absArg);
   }

   while (!queue.empty()) {
      RooAbsArg *node = queue.back();
      queue.pop_back();
      if (!visited.insert(node).second)
         continue;

      out->change(node, opMode);

      // Only follow value-client links: that is exactly the propagation path
      // used by RooAbsArg::setOperMode with mode==ADirty.
      if (opMode == RooAbsArg::ADirty) {
         for (auto *client : node->valueClients()) {
            queue.push_back(client);
         }
      }
   }
   return out;
}

void Evaluator::print(std::ostream &os)
{
   std::cout << "--- RooFit BatchMode evaluation ---\n";

   std::vector<int> widths{9, 37, 20, 9, 10, 20};

   auto printElement = [&](int iCol, auto const &t) {
      const char separator = ' ';
      os << separator << std::left << std::setw(widths[iCol]) << std::setfill(separator) << t;
      os << "|";
   };

   auto printHorizontalRow = [&]() {
      int n = 0;
      for (int w : widths) {
         n += w + 2;
      }
      for (int i = 0; i < n; i++) {
         os << '-';
      }
      os << "|\n";
   };

   printHorizontalRow();

   os << "|";
   printElement(0, "Index");
   printElement(1, "Name");
   printElement(2, "Class");
   printElement(3, "Size");
   printElement(4, "From Data");
   printElement(5, "1st value");
   std::cout << "\n";

   printHorizontalRow();

   for (std::size_t iNode = 0; iNode < _nodes.size(); ++iNode) {
      auto &nodeInfo = _nodes[iNode];
      RooAbsArg *node = nodeInfo.absArg;

      auto span = _evalContextCPU.at(node);

      os << "|";
      printElement(0, iNode);
      printElement(1, node->GetName());
      printElement(2, node->ClassName());
      printElement(3, nodeInfo.outputSize);
      printElement(4, nodeInfo.fromArrayInput);
      printElement(5, span[0]);

      std::cout << "\n";
   }

   printHorizontalRow();
}

/// Gets all the parameters of the RooAbsReal. This is in principle not
/// necessary, because we can always ask the RooAbsReal itself, but the
/// Evaluator has the cached information to get the answer quicker.
/// Therefore, this is not meant to be used in general, just where it matters.
/// \warning If we find another solution to get the parameters efficiently,
/// this function might be removed without notice.
RooArgSet Evaluator::getParameters() const
{
   RooArgSet parameters;
   for (auto &nodeInfo : _nodes) {
      if (nodeInfo.isValueServer && nodeInfo.absArg->isFundamental()) {
         parameters.add(*nodeInfo.absArg);
      }
   }
   // Just like in RooAbsArg::getParameters(), we sort the parameters alphabetically.
   parameters.sort();
   return parameters;
}

/// \brief Sets the offset mode for evaluation.
///
/// This function sets the offset mode for evaluation to the specified mode.
/// It updates the offset mode for both CPU and CUDA evaluation contexts.
///
/// \param mode The offset mode to be set.
///
/// \note This function marks reducer nodes as dirty if the offset mode is
///       changed, because only reducer nodes can use offsetting.
void Evaluator::setOffsetMode(RooFit::EvalContext::OffsetMode mode)
{
   if (mode == _evalContextCPU._offsetMode)
      return;

   _evalContextCPU._offsetMode = mode;
   _evalContextCUDA._offsetMode = mode;

   for (auto &nodeInfo : _nodes) {
      if (nodeInfo.absArg->isReducerNode()) {
         nodeInfo.isDirty = true;
      }
   }
}

} // namespace RooFit
