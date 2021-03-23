#ifndef RooFitDriver_h
#define RooFitDriver_h

#include "RooRealProxy.h"

#include "RunContext.h"

class RooAbsData;
class RooAbsReal;

class RooFitDriver : public RooAbsReal {

public:
   RooFitDriver(){};
   RooFitDriver(const char *name, const char *title, RooAbsReal &absReal, RooAbsData &data, RooArgSet const &normSet);
   RooFitDriver(const RooFitDriver &other, const char *name = 0);
   virtual TObject *clone(const char *newname) const override { return new RooFitDriver(*this, newname); }

   RooArgSet *getParameters(const RooArgSet *depList, Bool_t stripDisconnected = kTRUE) const override;

   virtual Double_t defaultErrorLevel() const override { return _absReal->defaultErrorLevel(); }

protected:
   RooAbsReal *_absReal = nullptr;
   RooAbsData *_data = nullptr;
   RooArgSet _normSet;

   double getValV(const RooArgSet *normalisationSet = nullptr) const override;

   double evaluate() const override;

   RooSpan<double> evaluateSpan(RooBatchCompute::RunContext &evalData, const RooArgSet *normSet) const override;

private:
   mutable RooBatchCompute::RunContext _runContext;

   ClassDefOverride(RooFitDriver, 1)
};

#endif
