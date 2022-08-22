#ifndef RooFit_RooFitCore_MyCudaAddition_h
#define RooFit_RooFitCore_MyCudaAddition_h

#include <RooAbsReal.h>
#include <RooListProxy.h>

class RooRealVar;
class RooArgList ;

class MyCudaAddition : public RooAbsReal {
public:

  MyCudaAddition(const char *name, const char *title, const RooArgList& sumSet) ;

  MyCudaAddition(const MyCudaAddition& other, const char* name = nullptr);
  TObject* clone(const char* newname) const override { return new MyCudaAddition(*this, newname); }

  void computeBatch(cudaStream_t*, double* output, size_t nEvents, RooFit::Detail::DataMap const&) const override;

protected:

  RooListProxy _set ;            ///< set of terms to be summed

  double evaluate() const override;

  inline bool canComputeBatchWithCuda() const override { return true; }
};

#endif
