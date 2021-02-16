#include "RooBatchCompute.h"

template<typename Function, typename ...Args_t>
RooSpan<double> callBatchCompute(RooAbsReal const* caller,
                                 Function function,
                                 RooBatchCompute::RunContext& runContext,
                                 Args_t const& ... arguments
                                 ) {
    return (RooBatchCompute::dispatch->*function) (caller, runContext, arguments...);
}
