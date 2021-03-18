// Author: Emmanouil Michalainas, CERN 3 March 2021

#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_CUH
#define ROOFIT_BATCHCOMPUTE_BATCHES_CUH

#include "RooSpan.h"

#include <stdint.h>
#include <vector>

class RooAbsReal;

namespace RooBatchCompute {

struct RunContext;

constexpr int nParams=8;

class Batch {
  private:
    double scalar=0;
    const double* __restrict array=nullptr;
    bool isVector=false;

  public:
    Batch() = default;
    Batch(const double* __restrict array, bool isVector)
      : scalar{array[0]}, array{array}, isVector{isVector}
    {}

    inline void set(double _scalar, const double* __restrict _array, bool _isVector)
    {
      scalar = _scalar;
      array = _array;
      isVector = _isVector;
    }
    __device__ constexpr double operator[](size_t i) const noexcept
    {
      return isVector ? array[i] : scalar;
    }    
}; //end class Batch


class Batches {
  private:
    size_t nEvents=0;
    Batch arrays[nParams];

  public:
    double* __restrict output=nullptr;

    Batches(const RooAbsReal* pdf, RunContext& evalData, const std::vector<const RooAbsReal*>& vars, const std::vector<RooSpan<const double>>& spans);
    __host__ __device__ constexpr size_t getNEvents()
    {
      return nEvents;
    }
    __device__ inline Batch operator[] (int batchIdx)
    {
      return arrays[batchIdx];
    }
}; //end class Batches

}

#endif /* ROOFIT_BATCHCOMPUTE_BATCHES_H */
