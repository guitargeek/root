#include "MyCudaAddition.h"

#include "Riostream.h"
#include "RooRealSumFunc.h"
#include "RooRealSumPdf.h"
#include "RooProduct.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNameReg.h"
#include "RooNLLVar.h"
#include "RooNLLVarNew.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
#include "RooBatchCompute.h"

#include <algorithm>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////
/// Constructor with a single set consisting of RooAbsReal.
/// \param[in] name Name of the PDF
/// \param[in] title Title
/// \param[in] sumSet The value of the function will be the sum of the values in this set

MyCudaAddition::MyCudaAddition(const char* name, const char* title, const RooArgList& sumSet)
  : RooAbsReal(name, title)
  , _set("!set","set of components",this)
{
  for (const auto comp : sumSet) {
    _set.add(*comp) ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

MyCudaAddition::MyCudaAddition(const MyCudaAddition& other, const char* name)
    : RooAbsReal(other, name)
    , _set("!set",this,other._set)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate and return current value of self

double MyCudaAddition::evaluate() const
{
  double sum(0);
  const RooArgSet* nset = _set.nset() ;

  for (auto* comp : static_range_cast<RooAbsReal*>(_set)) {
    const double tmp = comp->getVal(nset);
    sum += tmp ;
  }
  return sum ;
}



__global__
void additionKernel(std::size_t n, double const*x, double const*y, double *output)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i < n) output[i] = x[i] + y[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Compute addition of PDFs in batches.
void MyCudaAddition::computeBatch(cudaStream_t* stream, double* output, size_t nEvents, RooFit::Detail::DataMap const& dataMap) const
{
   auto xSpan = dataMap.at(&_set[0]);
   auto ySpan = dataMap.at(&_set[1]);

   if(stream) {
       // CUDA
       std::cout << "MyCudaAddition CUDA" << std::endl;
       int threadsPerBlock = 256;
       int blocksInGrid = std::ceil( double(nEvents) / threadsPerBlock );
       additionKernel<<<blocksInGrid, threadsPerBlock, 0, *stream>>>(nEvents, xSpan.data(), ySpan.data(), output);
   } else {
       // CPU
       std::cout << "MyCudaAddition CPU" << std::endl;
       for(std::size_t i = 0; i < nEvents; ++i) {
          output[i] = xSpan[i] + ySpan[i];
       }
   }
}
