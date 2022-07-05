#ifndef RooFit_HistFactory_GaussianConstraint_h
#define RooFit_HistFactory_GaussianConstraint_h

#include <RooAbsPdf.h>
#include <RooRealProxy.h>

//namespace RooStats {
//namespace HistFactory {

class GaussianConstraint : public RooAbsPdf {
public:
  GaussianConstraint() {}
  GaussianConstraint(const char *name, const char *title,
         RooAbsReal& x, RooAbsReal& mean, double sigma);
  GaussianConstraint(const GaussianConstraint& other, const char* name=nullptr);
  TObject* clone(const char* newname) const override {
    return new GaussianConstraint(*this,newname);
  }

  /// Get the x variable.
  RooAbsReal const& getX() const { return _x.arg(); }

  /// Get the mean parameter.
  RooAbsReal const& getMean() const { return _mean.arg(); }

  /// Get the sigma parameter.
  double getSigma() const { return _sigma; }

  /// Set the sigma parameter.
  double setSigma(double sigma) { return _sigma = sigma; }

  bool selfNormalized() const override { return true; }

private:

  RooRealProxy _x;
  RooRealProxy _mean;
  double _sigma;

  double evaluate() const override;

  ClassDefOverride(GaussianConstraint,1)
};

//}
//}

#endif
