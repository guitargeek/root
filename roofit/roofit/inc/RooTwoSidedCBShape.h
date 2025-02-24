#ifndef ROOT_Hfitter_RooTwoSidedCBShape
#define ROOT_Hfitter_RooTwoSidedCBShape

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"

class RooRealVar;


class RooTwoSidedCBShape : public RooAbsPdf {
public:
  RooTwoSidedCBShape() {}
  RooTwoSidedCBShape(const char *name, const char *title, RooAbsReal& _m,
               RooAbsReal& _m0, RooAbsReal& _sigma,
               RooAbsReal& _alphaLo, RooAbsReal& _nLo,
               RooAbsReal& _alphaHi, RooAbsReal& _nHi);

  RooTwoSidedCBShape(const RooTwoSidedCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooTwoSidedCBShape(*this,newname); }

  inline virtual ~RooTwoSidedCBShape() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  
protected:

  Double_t gaussianIntegral(Double_t tmin, Double_t tmax) const;
  Double_t powerLawIntegral(Double_t tmin, Double_t tmax, Double_t alpha, Double_t n) const;
  Double_t evaluate() const;

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaLo;
  RooRealProxy nLo;
  RooRealProxy alphaHi;
  RooRealProxy nHi;


private:

  ClassDef(RooTwoSidedCBShape,1)
};


#endif
