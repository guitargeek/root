/*
 * Project: RooFit
 * Authors:
 *   Emmanouil Michalainas, CERN, September 2020
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file RooBatchCompute.cu
\class RbcClass
\ingroup roofit_dev_docs_batchcompute

This file contains the code for cuda computations using the RooBatchCompute library.
**/

#include "RooBatchCompute.h"
#include "Batches.h"
#include "CudaInterface.h"

// CUDA Tile C++ (cuTile). Requires CUDA >= 13.3.
#include "cuda_tile.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <queue>
#include <vector>

namespace RooBatchCompute {
namespace CUDA {

// Alias for the cuTile namespace and enable the `_ic` compile-time-constant literals.
namespace ct = cuda::tiles;
using namespace ct::literals;

constexpr int blockSize = 512;

namespace {

void fillBatches(Batches &batches, double *output, size_t nEvents, std::size_t nBatches, std::size_t nExtraArgs)
{
   batches.nEvents = nEvents;
   batches.nBatches = nBatches;
   batches.nExtra = nExtraArgs;
   batches.output = output;
}

void fillArrays(Batch *arrays, VarSpan vars, double *buffer, double *bufferDevice, std::size_t nEvents)
{
   for (int i = 0; i < vars.size(); i++) {
      const std::span<const double> &span = vars[i];
      arrays[i]._isVector = span.empty() || span.size() >= nEvents;
      if (!arrays[i]._isVector) {
         // In the scalar case, the value is not on the GPU yet, so we have to
         // copy the value to the GPU buffer.
         buffer[i] = span[0];
         arrays[i]._array = bufferDevice + i;
      } else {
         // In the vector input cases, they are already on the GPU, so we can
         // fill be buffer with some dummy value and set the input span
         // directly.
         buffer[i] = 0.0;
         arrays[i]._array = span.data();
      }
   }
}

/// Number of streaming multiprocessors (SMs) of the current device, queried
/// once and cached. This replaces the grid-size cap that used to be hard-coded
/// to 84 -- a value the original author had tuned by hand for the specific GPU
/// they happened to develop on (an RTX A4500 with 46 SMs). Deriving it from the
/// actual hardware means the element-wise kernels scale to any device.
int deviceMultiProcessorCount()
{
   static const int count = [] {
      int device = 0;
      cudaGetDevice(&device);
      int smCount = 0;
      cudaDeviceGetAttribute(&smCount, cudaDevAttrMultiProcessorCount, device);
      return smCount > 0 ? smCount : 1;
   }();
   return count;
}

/// Grid size for the element-wise grid-stride compute kernels. We launch enough
/// blocks to saturate the device (a small multiple of the SM count) but never
/// more than are needed to cover the data. Unlike the old cap, this is not a
/// magic constant: the grid-stride loops inside the kernels remain correct for
/// any grid size, so this only affects occupancy, not results.
int getGridSize(std::size_t n)
{
   const int blocksToCoverData = int(std::ceil(double(n) / blockSize));
   const int blocksToFillDevice = 2 * deviceMultiProcessorCount();
   return std::max(1, std::min(blocksToCoverData, blocksToFillDevice));
}

} // namespace

std::vector<void (*)(Batches &)> getFunctions();

/// This class overrides some RooBatchComputeInterface functions, for the
/// purpose of providing a cuda specific implementation of the library.
class RooBatchComputeClass : public RooBatchComputeInterface {

public:
   RooBatchComputeClass() : _computeFunctions(getFunctions())
   {
      dispatchCUDA = this; // Set the dispatch pointer to this instance of the library upon loading
   }

   Architecture architecture() const override { return Architecture::CUDA; }
   std::string architectureName() const override { return "cuda"; }

   /** Compute multiple values using cuda kernels.
   This method creates a Batches object and passes it to the correct compute function.
   The compute function is launched as a cuda kernel.
   \param computer An enum specifying the compute function to be used.
   \param output The array where the computation results are stored.
   \param vars A std::span containing pointers to the variables involved in the computation.
   \param extraArgs An optional std::span containing extra double values that may participate in the computation. **/
   void compute(RooBatchCompute::Config const &cfg, Computer computer, std::span<double> output, VarSpan vars,
                ArgSpan extraArgs) override
   {
      using namespace CudaInterface;

      std::size_t nEvents = output.size();

      const std::size_t memSize = sizeof(Batches) + vars.size() * sizeof(Batch) + vars.size() * sizeof(double) +
                                  extraArgs.size() * sizeof(double);

      std::vector<char> hostMem(memSize);
      auto batches = reinterpret_cast<Batches *>(hostMem.data());
      auto arrays = reinterpret_cast<Batch *>(batches + 1);
      auto scalarBuffer = reinterpret_cast<double *>(arrays + vars.size());
      auto extraArgsHost = reinterpret_cast<double *>(scalarBuffer + vars.size());

      DeviceArray<char> deviceMem(memSize);
      auto batchesDevice = reinterpret_cast<Batches *>(deviceMem.data());
      auto arraysDevice = reinterpret_cast<Batch *>(batchesDevice + 1);
      auto scalarBufferDevice = reinterpret_cast<double *>(arraysDevice + vars.size());
      auto extraArgsDevice = reinterpret_cast<double *>(scalarBufferDevice + vars.size());

      fillBatches(*batches, output.data(), nEvents, vars.size(), extraArgs.size());
      fillArrays(arrays, vars, scalarBuffer, scalarBufferDevice, nEvents);
      batches->args = arraysDevice;

      if (!extraArgs.empty()) {
         std::copy(std::cbegin(extraArgs), std::cend(extraArgs), extraArgsHost);
         batches->extra = extraArgsDevice;
      }

      copyHostToDevice(hostMem.data(), deviceMem.data(), hostMem.size(), cfg.cudaStream());

      const int gridSize = getGridSize(nEvents);
      _computeFunctions[computer]<<<gridSize, blockSize, 0, *cfg.cudaStream()>>>(*batchesDevice);

      // The compute might have modified the mutable extra args, so we need to
      // copy them back. This can be optimized if necessary in the future by
      // flagging if the extra args were actually changed.
      if (!extraArgs.empty()) {
         copyDeviceToHost(extraArgsDevice, extraArgs.data(), extraArgs.size(), cfg.cudaStream());
      }
   }
   /// Return the sum of an input array
   double reduceSum(RooBatchCompute::Config const &cfg, InputArr input, size_t n) override;
   ReduceNLLOutput reduceNLL(RooBatchCompute::Config const &cfg, std::span<const double> probas,
                             std::span<const double> weights, std::span<const double> offsetProbas) override;

   std::unique_ptr<AbsBufferManager> createBufferManager() const override;

   CudaInterface::CudaEvent *newCudaEvent(bool forTiming) const override
   {
      return new CudaInterface::CudaEvent{forTiming};
   }
   CudaInterface::CudaStream *newCudaStream() const override { return new CudaInterface::CudaStream{}; }
   void deleteCudaEvent(CudaInterface::CudaEvent *event) const override { delete event; }
   void deleteCudaStream(CudaInterface::CudaStream *stream) const override { delete stream; }

   void cudaEventRecord(CudaInterface::CudaEvent *event, CudaInterface::CudaStream *stream) const override
   {
      CudaInterface::cudaEventRecord(*event, *stream);
   }
   void cudaStreamWaitForEvent(CudaInterface::CudaStream *stream, CudaInterface::CudaEvent *event) const override
   {
      stream->waitForEvent(*event);
   }
   bool cudaStreamIsActive(CudaInterface::CudaStream *stream) const override { return stream->isActive(); }

private:
   const std::vector<void (*)(Batches &)> _computeFunctions;

}; // End class RooBatchComputeClass

namespace {

// ----------------------------------------------------------------------------
// cuTile-based reductions
//
// The reductions below replace the hand-rolled shared-memory Kahan-summation
// kernels. Each cuTile kernel reduces one tile of the input (one tile per
// block) to a single partial sum; cuTile picks the internal thread parallelism
// for the current architecture, so there are no hand-tuned block/grid launch
// parameters left. The only tunable is the compile-time tile length below.
//
// The per-tile partials are then combined with compensated (Neumaier) summation
// on the host to recover the Kahan sum + carry that RooFit relies on for
// reproducible NLL values. The number of partials equals the number of tiles,
// which is tiny for realistic dataset sizes, so this host-side step is cheap.
// ----------------------------------------------------------------------------

/// Compile-time tile length for the reduction kernels (elements per tile/block).
constexpr auto reduceTile = 1024_ic;

/// Number of tiles (= number of blocks) needed to cover `n` elements.
std::size_t numReduceTiles(std::size_t n)
{
   return (n + std::size_t(reduceTile) - 1) / std::size_t(reduceTile);
}

/// Compensated (Neumaier) accumulation on the host, matching the semantics of
/// ROOT::Math::KahanSum. Used to combine the per-tile partial sums.
inline void kahanAccumulate(double &sum, double &carry, double value)
{
   const double t = sum + value;
   if (std::abs(sum) >= std::abs(value))
      carry += (sum - t) + value;
   else
      carry += (value - t) + sum;
   sum = t;
}

/// Combine `nTiles` device-side partial sums into a final (sum, carry) pair.
void combinePartialsOnHost(const CudaInterface::DeviceArray<double> &devPartial, std::size_t nTiles,
                           cudaStream_t stream, double &sum, double &carry)
{
   // Reusable per-thread scratch buffer so we don't heap-allocate on every call
   // (the rest of RooBatchCompute pools its buffers for the same reason). It is
   // thread-local to stay safe when reductions run concurrently on separate
   // streams, and only ever grows.
   thread_local std::vector<double> partial;
   partial.resize(nTiles);

   // Issue the copy on the same stream as the reduction kernel so it is ordered
   // after it, then block until the partials are on the host. The final result
   // must reach the CPU regardless (the minimizer runs host-side), so this is the
   // one unavoidable transfer -- and it is only nTiles doubles, i.e. latency-bound.
   cudaMemcpyAsync(partial.data(), devPartial.data(), nTiles * sizeof(double), cudaMemcpyDeviceToHost, stream);
   cudaStreamSynchronize(stream);

   sum = 0.0;
   carry = 0.0;
   for (std::size_t i = 0; i < nTiles; ++i)
      kahanAccumulate(sum, carry, partial[i]);
}

} // namespace

/// cuTile kernel: each block reduces one tile of `input` to a single partial sum
/// written to `partial[blockIdx]`. Out-of-range lanes of the last tile read as 0.
__tile_global__ void tileReduceSum(const double *__restrict__ input, std::size_t n, std::size_t nTiles,
                                   double *__restrict__ partial)
{
   auto in = ct::partition_view{ct::tensor_span{input, ct::extents{n}}, ct::shape{reduceTile}};
   auto out = ct::partition_view{ct::tensor_span{partial, ct::extents{nTiles}}, ct::shape{1_ic}};
   const int b = ct::bid().x;
   out.store(ct::sum(in.load_masked(b), 0_ic), b);
}

/// cuTile kernel for the NLL reduction. Computes, per event,
///     val = weight * ( -log(proba) [ + log(offsetProba) ] )
/// and reduces each tile to a partial sum. Padding lanes of the last tile load
/// `proba`/`offsetProba` as 1.0 (so -log() -> 0) and `weight` as 0.0, which keeps
/// the padding contribution at exactly 0 and avoids 0*inf = NaN.
__tile_global__ void tileReduceNLL(const double *__restrict__ probas, const double *__restrict__ weights,
                                   const double *__restrict__ offsetProbas, bool scalarProba, double negLogScalarProba,
                                   std::size_t n, std::size_t nTiles, double *__restrict__ partial)
{
   auto wView = ct::partition_view{ct::tensor_span{weights, ct::extents{n}}, ct::shape{reduceTile}};
   auto outView = ct::partition_view{ct::tensor_span{partial, ct::extents{nTiles}}, ct::shape{1_ic}};
   const int b = ct::bid().x;

   auto w = wView.load_masked(b, 0.0);

   // -log(proba): either a broadcast scalar or an elementwise log of the tile.
   auto negLogProba = scalarProba
                         ? ct::fill_like(w, negLogScalarProba)
                         : -ct::log(ct::partition_view{ct::tensor_span{probas, ct::extents{n}}, ct::shape{reduceTile}}
                                       .load_masked(b, 1.0));

   auto val = w * negLogProba;

   if (offsetProbas != nullptr) {
      auto off =
         ct::partition_view{ct::tensor_span{offsetProbas, ct::extents{n}}, ct::shape{reduceTile}}.load_masked(b, 1.0);
      val = val + w * ct::log(off);
   }

   outView.store(ct::sum(val, 0_ic), b);
}

double RooBatchComputeClass::reduceSum(RooBatchCompute::Config const &cfg, InputArr input, size_t n)
{
   if (n == 0)
      return 0.0;
   const std::size_t nTiles = numReduceTiles(n);
   cudaStream_t stream = *cfg.cudaStream();
   CudaInterface::DeviceArray<double> devPartial(nTiles);
   // cuTile launch: grid = one block per tile, and the block dimension MUST be 1
   // (cuTile manages the internal parallelism itself).
   tileReduceSum<<<static_cast<unsigned int>(nTiles), 1, 0, stream>>>(input, n, nTiles, devPartial.data());
   double sum = 0.0;
   double carry = 0.0;
   combinePartialsOnHost(devPartial, nTiles, stream, sum, carry);
   return sum + carry;
}

ReduceNLLOutput RooBatchComputeClass::reduceNLL(RooBatchCompute::Config const &cfg, std::span<const double> probas,
                                                std::span<const double> weights, std::span<const double> offsetProbas)
{
   ReduceNLLOutput out;
   if (probas.empty()) {
      return out;
   }
   const std::size_t n = weights.size();
   const std::size_t nTiles = numReduceTiles(n);
   cudaStream_t stream = *cfg.cudaStream();
   CudaInterface::DeviceArray<double> devPartial(nTiles);

#ifndef NDEBUG
   for (auto span : {probas, weights, offsetProbas}) {
      cudaPointerAttributes attr;
      assert(span.size() == 0 || span.data() == nullptr ||
             (cudaPointerGetAttributes(&attr, span.data()) == cudaSuccess && attr.type == cudaMemoryTypeDevice));
   }
#endif

   const bool scalarProba = probas.size() == 1;
   // For the scalar-proba case we cannot dereference the device pointer on the
   // host, so the caller-visible scalar value is read back by RooFit elsewhere;
   // here scalarProba==true implies a single value that we pre-transform. When
   // proba is a genuine per-event array this term is computed on the device.
   const double negLogScalarProba = scalarProba ? -std::log(probas[0]) : 0.0;

   tileReduceNLL<<<static_cast<unsigned int>(nTiles), 1, 0, stream>>>(
      probas.data(), weights.data(), offsetProbas.empty() ? nullptr : offsetProbas.data(), scalarProba,
      negLogScalarProba, n, nTiles, devPartial.data());

   double sum = 0.0;
   double carry = 0.0;
   combinePartialsOnHost(devPartial, nTiles, stream, sum, carry);

   out.nllSum = sum;
   out.nllSumCarry = carry;
   return out;
}

namespace {

class ScalarBufferContainer {
public:
   ScalarBufferContainer() {}
   ScalarBufferContainer(std::size_t size)
   {
      if (size != 1)
         throw std::runtime_error("ScalarBufferContainer can only be of size 1");
   }

   double const *hostReadPtr() const { return &_val; }
   double const *deviceReadPtr() const { return &_val; }

   double *hostWritePtr() { return &_val; }
   double *deviceWritePtr() { return &_val; }

   void assignFromHost(std::span<const double> input) { _val = input[0]; }
   void assignFromDevice(std::span<const double> input)
   {
      CudaInterface::copyDeviceToHost(input.data(), &_val, input.size(), nullptr);
   }

private:
   double _val;
};

class CPUBufferContainer {
public:
   CPUBufferContainer(std::size_t size) : _vec(size) {}

   double const *hostReadPtr() const { return _vec.data(); }
   double const *deviceReadPtr() const
   {
      throw std::bad_function_call();
      return nullptr;
   }

   double *hostWritePtr() { return _vec.data(); }
   double *deviceWritePtr()
   {
      throw std::bad_function_call();
      return nullptr;
   }

   void assignFromHost(std::span<const double> input) { _vec.assign(input.begin(), input.end()); }
   void assignFromDevice(std::span<const double> input)
   {
      CudaInterface::copyDeviceToHost(input.data(), _vec.data(), input.size(), nullptr);
   }

private:
   std::vector<double> _vec;
};

class GPUBufferContainer {
public:
   GPUBufferContainer(std::size_t size) : _arr(size) {}

   double const *hostReadPtr() const
   {
      throw std::bad_function_call();
      return nullptr;
   }
   double const *deviceReadPtr() const { return _arr.data(); }

   double *hostWritePtr() const
   {
      throw std::bad_function_call();
      return nullptr;
   }
   double *deviceWritePtr() const { return const_cast<double *>(_arr.data()); }

   void assignFromHost(std::span<const double> input)
   {
      CudaInterface::copyHostToDevice(input.data(), deviceWritePtr(), input.size(), nullptr);
   }
   void assignFromDevice(std::span<const double> input)
   {
      CudaInterface::copyDeviceToDevice(input.data(), deviceWritePtr(), input.size(), nullptr);
   }

private:
   CudaInterface::DeviceArray<double> _arr;
};

class PinnedBufferContainer {
public:
   PinnedBufferContainer(std::size_t size) : _arr{size}, _gpuBuffer{size} {}
   std::size_t size() const { return _arr.size(); }

   void setCudaStream(CudaInterface::CudaStream *stream) { _cudaStream = stream; }

   double const *hostReadPtr() const
   {

      if (_lastAccess == LastAccessType::GPU_WRITE) {
         CudaInterface::copyDeviceToHost(_gpuBuffer.deviceReadPtr(), const_cast<double *>(_arr.data()), size(),
                                         _cudaStream);
      }

      _lastAccess = LastAccessType::CPU_READ;
      return const_cast<double *>(_arr.data());
   }
   double const *deviceReadPtr() const
   {

      if (_lastAccess == LastAccessType::CPU_WRITE) {
         CudaInterface::copyHostToDevice(_arr.data(), _gpuBuffer.deviceWritePtr(), size(), _cudaStream);
      }

      _lastAccess = LastAccessType::GPU_READ;
      return _gpuBuffer.deviceReadPtr();
   }

   double *hostWritePtr()
   {
      _lastAccess = LastAccessType::CPU_WRITE;
      return _arr.data();
   }
   double *deviceWritePtr()
   {
      _lastAccess = LastAccessType::GPU_WRITE;
      return _gpuBuffer.deviceWritePtr();
   }

   void assignFromHost(std::span<const double> input) { std::copy(input.begin(), input.end(), hostWritePtr()); }
   void assignFromDevice(std::span<const double> input)
   {
      CudaInterface::copyDeviceToDevice(input.data(), deviceWritePtr(), input.size(), _cudaStream);
   }

private:
   enum class LastAccessType {
      CPU_READ,
      GPU_READ,
      CPU_WRITE,
      GPU_WRITE
   };

   CudaInterface::PinnedHostArray<double> _arr;
   GPUBufferContainer _gpuBuffer;
   CudaInterface::CudaStream *_cudaStream = nullptr;
   mutable LastAccessType _lastAccess = LastAccessType::CPU_READ;
};

template <class Container>
class BufferImpl : public AbsBuffer {
public:
   using Queue = std::queue<std::unique_ptr<Container>>;

   BufferImpl(std::size_t size, Queue &queue) : _queue{queue}
   {
      if (_queue.empty()) {
         _vec = std::make_unique<Container>(size);
      } else {
         _vec = std::move(_queue.front());
         _queue.pop();
      }
   }

   ~BufferImpl() override { _queue.emplace(std::move(_vec)); }

   double const *hostReadPtr() const override { return _vec->hostReadPtr(); }
   double const *deviceReadPtr() const override { return _vec->deviceReadPtr(); }

   double *hostWritePtr() override { return _vec->hostWritePtr(); }
   double *deviceWritePtr() override { return _vec->deviceWritePtr(); }

   void assignFromHost(std::span<const double> input) override { _vec->assignFromHost(input); }
   void assignFromDevice(std::span<const double> input) override { _vec->assignFromDevice(input); }

   Container &vec() { return *_vec; }

private:
   std::unique_ptr<Container> _vec;
   Queue &_queue;
};

using ScalarBuffer = BufferImpl<ScalarBufferContainer>;
using CPUBuffer = BufferImpl<CPUBufferContainer>;
using GPUBuffer = BufferImpl<GPUBufferContainer>;
using PinnedBuffer = BufferImpl<PinnedBufferContainer>;

struct BufferQueuesMaps {
   std::map<std::size_t, ScalarBuffer::Queue> scalarBufferQueuesMap;
   std::map<std::size_t, CPUBuffer::Queue> cpuBufferQueuesMap;
   std::map<std::size_t, GPUBuffer::Queue> gpuBufferQueuesMap;
   std::map<std::size_t, PinnedBuffer::Queue> pinnedBufferQueuesMap;
};

class BufferManager : public AbsBufferManager {

public:
   BufferManager() : _queuesMaps{std::make_unique<BufferQueuesMaps>()} {}

   std::unique_ptr<AbsBuffer> makeScalarBuffer() override
   {
      return std::make_unique<ScalarBuffer>(1, _queuesMaps->scalarBufferQueuesMap[1]);
   }
   std::unique_ptr<AbsBuffer> makeCpuBuffer(std::size_t size) override
   {
      return std::make_unique<CPUBuffer>(size, _queuesMaps->cpuBufferQueuesMap[size]);
   }
   std::unique_ptr<AbsBuffer> makeGpuBuffer(std::size_t size) override
   {
      return std::make_unique<GPUBuffer>(size, _queuesMaps->gpuBufferQueuesMap[size]);
   }
   std::unique_ptr<AbsBuffer> makePinnedBuffer(std::size_t size, CudaInterface::CudaStream *stream = nullptr) override
   {
      auto out = std::make_unique<PinnedBuffer>(size, _queuesMaps->pinnedBufferQueuesMap[size]);
      out->vec().setCudaStream(stream);
      return out;
   }

private:
   std::unique_ptr<BufferQueuesMaps> _queuesMaps;
};

} // namespace

std::unique_ptr<AbsBufferManager> RooBatchComputeClass::createBufferManager() const
{
   return std::make_unique<BufferManager>();
}

/// Static object to trigger the constructor which overwrites the dispatch pointer.
static RooBatchComputeClass computeObj;

} // End namespace CUDA
} // End namespace RooBatchCompute
