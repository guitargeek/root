#ifndef RooFit_RooCruijff_h
#define RooFit_RooCruijff_h

#include <RooAbsPdf.h>
#include <RooConstVar.h>
#include <RooRealProxy.h>

class RooCruijff : public RooAbsPdf {
public:
   RooCruijff() {}
   inline RooCruijff(const char *name, const char *title, RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigmaL,
                     RooAbsReal &sigmaR, RooAbsReal &alphaL, RooAbsReal &alphaR)
      : RooCruijff{name, title, x, mu, sigmaL, sigmaR, alphaL, alphaR, RooFit::RooConst(0.0)}
   {
   }
   RooCruijff(const char *name, const char *title, RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigmaL,
              RooAbsReal &sigmaR, RooAbsReal &alphaL, RooAbsReal &alphaR, RooAbsReal &beta);
   RooCruijff(const RooCruijff &other, const char *name = nullptr);
   TObject *clone(const char *newname) const override { return new RooCruijff{*this, newname}; }

protected:
   RooRealProxy _x;
   RooRealProxy _mu;
   RooRealProxy _sigmaL;
   RooRealProxy _sigmaR;
   RooRealProxy _alphaL;
   RooRealProxy _alphaR;
   RooRealProxy _beta;

   double evaluate() const override;
   void computeBatch(double *output, size_t size, RooFit::Detail::DataMap const &) const override;
   inline bool canComputeBatchWithCuda() const override { return true; }

private:
   ClassDefOverride(RooCruijff, 1); // Cruijff PDF
};

#endif
