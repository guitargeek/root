// Author: Emmanouil Michalainas, CERN 3 March 2021

#include "Batches.h"
#include "RunContext.h"

#include <algorithm>

RooBatchCompute::Batches::Batches(const RooAbsReal* pdf, RunContext& evalData, const std::vector<RooSpan<const double>>& spans, double stackArr[nParams][bufferSize])
{
  for (size_t i=0; i<spans.size(); i++)
    if (spans[i].size()>1)
    {
      nEvents = spans[i].size();
      arrays[i] = spans[i].data();
    }
    else
    {
      std::fill_n(stackArr[i], bufferSize, spans[i][0]);
      arrays[i] = stackArr[i];
    }
  output = evalData.makeBatch(pdf, nEvents).data();
}
