#include <RooStats/HistFactory/GaussianConstraint.h>

#include <RooMath.h>

ClassImp(GaussianConstraint);


GaussianConstraint::GaussianConstraint(const char *name, const char *title,
          RooAbsReal& x, RooAbsReal& mean,
          double sigma) :
  RooAbsPdf{name,title},
  _x{"x","Observable",this,x},
  _mean{"mean","Mean",this,mean},
  _sigma{sigma}
{
}


GaussianConstraint::GaussianConstraint(const GaussianConstraint& other, const char* name) :
  RooAbsPdf{other,name}, _x{"x",this,other._x}, _mean{"mean",this,other._mean},
  _sigma{other._sigma}
{
}


double GaussianConstraint::evaluate() const
{
  constexpr double oneOverSqrtTwoPi = 1./std::sqrt(TMath::TwoPi());
  const double arg = _x - _mean;
  return oneOverSqrtTwoPi / _sigma * std::exp(-0.5*arg*arg/(_sigma*_sigma));
}
