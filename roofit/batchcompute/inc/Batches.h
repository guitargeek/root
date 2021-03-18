// Author: Emmanouil Michalainas, CERN 3 March 2021

#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_H
#define ROOFIT_BATCHCOMPUTE_BATCHES_H

#include "RooSpan.h"

#include <stdint.h>
#include <vector>

class RooAbsReal;

namespace RooBatchCompute {

struct RunContext;

constexpr int nParams=8;
constexpr int bufferSize=64;

class Batches {
  private:
    size_t nEvents=0;
    const double* __restrict arrays[nParams] {nullptr};

  public:
    double* __restrict output=nullptr;
    
    Batches(const RooAbsReal* pdf, RunContext& evalData, const std::vector<RooSpan<const double>>& spans, double stackArr[nParams][bufferSize]);
    inline size_t setNEvents(size_t n=bufferSize)
    {
      size_t temp = nEvents;
      nEvents = n;
      return temp;
    }
    constexpr size_t getNEvents()
    {
      return nEvents;
    }
    constexpr const double* operator[] (int batchIdx)
    {
      return arrays[batchIdx];
    }
    inline void advance(const std::vector<RooSpan<const double>>& vars)
    {
      for (size_t i=0; i<vars.size(); i++)
        arrays[i] += (vars[i].size()>1)*bufferSize;
      output += bufferSize;
    }
}; //end class Batches
} //end namespace RooBatchCompute

#endif // #ifdef ROOFIT_BATCHCOMPUTE_BATCHES_H 
