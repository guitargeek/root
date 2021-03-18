// RooBatchCompute library created September 2020 by Emmanouil Michalainas

#include "RooBatchCompute.h"
#include "Batches.h"

#define ROOBATCHCOMPUTE_BEGIN  0
#define ROOBATCHCOMPUTE_STEP   1

namespace RooBatchCompute {

class RooBatchComputeClass : public RooBatchComputeInterface {

  public:
    RooBatchComputeClass()
    {
      // Set the dispatch pointer to this instance of the library upon loading
      RooBatchCompute::dispatch = this;
    }

    void compute(const RooAbsReal* caller, Computer computer, RunContext& evalData, const std::vector<const RooAbsReal*>& , const std::vector<RooSpan<const double>>& spans) override  
    {
      double buffer[nParams][bufferSize];
      Batches batches(caller, evalData, spans, buffer);
      size_t nEvents = batches.setNEvents();
      while (nEvents > bufferSize) {
        computeFunctions[computer](batches);
        batches.advance(spans);
        nEvents -= bufferSize;
      }
      batches.setNEvents(nEvents);
      computeFunctions[computer](batches);
    }
    
}; // End class RooBatchComputeClass

/// Static object to trigger the constructor which overwrites the dispatch pointer.
static RooBatchComputeClass computeObj;

} //End namespace RooBatchCompute
