// Author: Emmanouil Michalainas, CERN 3 March 2021

#include "Batches.cuh"
#include "RunContext.h"

#include <iostream>
#include <new>

namespace {

double* errorCheckingCudaMalloc(size_t size)
{
  double* ret = nullptr;
  cudaError_t error = cudaMalloc(&ret, size);
  if (error != cudaSuccess)
  {
    std::cerr << "In " << __func__ << "(), " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(error) << std::endl;
    throw std::bad_alloc();
  }
  else return ret;
}

} // end namespace

namespace RooBatchCompute {

Batches::Batches(const RooAbsReal* pdf, RunContext& evalData, const std::vector<const RooAbsReal*>& vars, const std::vector<RooSpan<const double>>& spans)
{  
  for (int i=0; i<vars.size(); i++) {
    auto size = spans[i].size();
    if (size==1) {
      arrays[i].set(spans[i][0], nullptr, false);
      continue;
    }
    nEvents = size;
    // convert to bytes for allocation & copy
    size *= sizeof(double);
    
    double*& pointer = evalData.ownedMemoryCuda[vars[i]];
    if (pointer == nullptr) {
      pointer = errorCheckingCudaMalloc(size);
    }
    if (evalData.spansCuda[vars[i]] == nullptr) {
      evalData.spansCuda[vars[i]] = pointer;
      cudaMemcpy(pointer, spans[i].data(), size, cudaMemcpyHostToDevice);
    }
    arrays[i].set(0.0, pointer, true);
  }

  double*& pointer = evalData.ownedMemoryCuda[pdf];
  if (pointer == nullptr) {
    pointer = errorCheckingCudaMalloc(nEvents*sizeof(double));
  }
  evalData.spansCuda[pdf] = output = pointer;
}

} // end namespace RooBatchCompute
