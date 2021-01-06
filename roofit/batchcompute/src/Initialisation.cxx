#include "RooBatchCompute.h"

#include "TEnv.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <exception>


// First initialisation of the pointer. When implementations of the batch compute library are loaded,
// they will overwrite the pointer.
RooBatchCompute::RooBatchComputeInterface* RooBatchCompute::dispatch=nullptr;

namespace {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Dynamically load a library and throw exception in case of failure
void loadWithErrorChecking(std::string libName)
{
  const auto returnValue = gSystem->Load(libName.c_str());
  if (returnValue == -1 || returnValue == -2) {
    throw std::runtime_error("RooFit was unable to load its computation library " + libName);
  } else if (returnValue == 1) {
    // Library should not have been loaded before we tried to do it.
    throw std::logic_error("RooFit computation library " + libName + " was loaded before RooFit initialisation began.");
  } 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inspect hardware capabilities, and load the optimal library for RooFit computations.
void loadComputeLibrary() {

  // Check if user has disabled optimized libraries in .rootrc 
  if (gEnv->GetValue("RooFit.LoadOptimisedComputationLibrary",1) == 0) {
    loadWithErrorChecking("libRooBatchCompute_GENERIC");
    return;
  }

  // Check if cuda is enabled and available
  if (gEnv->GetValue("RooFit.Cuda",1)==1 && gSystem->Load("libcuda")>=0 && gSystem->Load("libRooBatchCompute_CUDA")==0 && RooBatchCompute::dispatch!=nullptr) {
    return;
  } else if (gDebug>0) {
    std::cout << "In " << __func__ << "(), " << __FILE__ << ":" << __LINE__ << ": Failed to load cuda implementation, trying cpu optimised implementations." << std::endl;
  }

#ifdef R__RF_ARCHITECTURE_SPECIFIC_LIBS
  
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx512cd") && __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512dq"))  {
    loadWithErrorChecking("libRooBatchCompute_AVX512");
    return;
  } else if (__builtin_cpu_supports("avx2")) {
    loadWithErrorChecking("libRooBatchCompute_AVX2");
    return;
  } else if (__builtin_cpu_supports("avx")) {
    loadWithErrorChecking("libRooBatchCompute_AVX");
    return;
  } else if (__builtin_cpu_supports("sse4.1")) {
    loadWithErrorChecking("libRooBatchCompute_SSE4.1");
    return;
  }
  
#endif //R__RF_ARCHITECTURE_SPECIFIC_LIBS

  if (gDebug>0) {
    std::cout << "In " << __func__ << "(), " << __FILE__ << ":" << __LINE__ << ": Vector instruction sets not supported, using generic implementation." << std::endl;
  }
  loadWithErrorChecking("libRooBatchCompute_GENERIC");

}

} //end anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// A RAII that performs RooFit's static initialisation.
static struct RooBatchComputeInitialiser {
  RooBatchComputeInitialiser() {
    loadComputeLibrary();
  }
} __RooBatchComputeInitialiser;

