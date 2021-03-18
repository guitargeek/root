// RooBatchCompute library created September 2020 by Emmanouil Michalainas

#include "RooBatchCompute.h"
#include "Batches.cuh"

#include "TEnv.h"

#include <iostream>

namespace RooBatchCompute {

class RooBatchComputeClass : public RooBatchComputeInterface {

  public:
    RooBatchComputeClass()
    {
      int nDevices=0;
      cudaError_t err = cudaGetDeviceCount(&nDevices);
      if (err!=cudaSuccess) {
        if (gDebug>0) {
          std::cerr << "In " << __func__ << "(), " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl;
        }
      } else if (nDevices>0)
      // Set the dispatch pointer to this instance of the library upon loading
      RooBatchCompute::dispatch = this;
    }

    void compute(const RooAbsReal* caller, Computer computer, RunContext& evalData, const std::vector<const RooAbsReal*>& vars, const std::vector<RooSpan<const double>>& spans) override  
    {
      Batches batches(caller, evalData, vars, spans);
      computeFunctions[computer]<<<1,1>>>(batches);
    }
    
    
}; // End class RooBatchComputeClass

/// Static object to trigger the constructor which overwrites the dispatch pointer.
static RooBatchComputeClass computeObj;

} //End namespace RooBatchCompute
