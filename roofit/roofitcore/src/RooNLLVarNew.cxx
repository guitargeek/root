#include "RooNLLVarNew.h"
#include "RooBatchCompute.h"

#include <numeric>
#include <stdexcept>
#include <vector>

ClassImp(RooNLLVarNew);

RooNLLVarNew::RooNLLVarNew(const char *name, const char *title, RooAbsPdf &pdf)
   : RooAbsReal(name, title), _pdf{"pdf", "pdf", this, pdf}
{
}

RooNLLVarNew::RooNLLVarNew(const RooNLLVarNew &other, const char *name)
   : RooAbsReal(other, name), _pdf{"pdf", this, other._pdf}
{
}

double RooNLLVarNew::getValV(const RooArgSet *) const
{
   throw std::runtime_error("RooNLLVarNew::getValV was called directly which should not happen!");
}

double RooNLLVarNew::evaluate() const
{
   throw std::runtime_error("RooNLLVarNew::evaluate was called directly which should not happen!");
}

RooSpan<double> RooNLLVarNew::evaluateSpan(RooBatchCompute::RunContext &, const RooArgSet *) const
{
   throw std::runtime_error("RooNLLVarNew::evaluatSpan was called directly which should not happen!");
}

RooSpan<const double> RooNLLVarNew::getValues(RooBatchCompute::RunContext &evalData, const RooArgSet *normSet) const
{

   if (!normSet) {
      throw std::runtime_error("RooNLLVarNew::getValues called without normalization set!");
   }

   auto item = evalData.spans.find(this);
   if (item != evalData.spans.end()) {
      return item->second;
   }

   auto values = _pdf->getValues(evalData, normSet);
   std::vector<double> logValues;
   logValues.resize(values.size());

   for (std::size_t i = 0; i < values.size(); ++i) { // CHECK_VECTORISE
      logValues[i] = RooBatchCompute::fast_log(values[i]);
   }

   auto outputData = evalData.makeBatch(this, 1);
   outputData[0] = -std::accumulate(logValues.begin(), logValues.end(), 0.0);

   return outputData;
}

RooArgSet *RooNLLVarNew::getParameters(const RooArgSet *depList, Bool_t stripDisconnected) const
{
   return _pdf->getParameters(depList, stripDisconnected);
}
