#include <RooCruijff.h>

#include <RooBatchCompute.h>

#include <RooFit/Detail/EvaluateFuncs.h>

ClassImp(RooCruijff);

RooCruijff::RooCruijff(const char *name, const char *title, RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigmaL,
                       RooAbsReal &sigmaR, RooAbsReal &alphaL, RooAbsReal &alphaR, RooAbsReal &beta)
   : RooAbsPdf{name, title},
     _x{"x", "Observable", this, x},
     _mu{"mu", "Mean", this, mu},
     _sigmaL{"sigmaL", "Width left", this, sigmaL},
     _sigmaR{"sigmaR", "Width right", this, sigmaR},
     _alphaL{"alphaL", "Shape left", this, alphaL},
     _alphaR{"alphaR", "Shape right", this, alphaR},
     _beta{"beta", "Beta", this, beta}
{
}

RooCruijff::RooCruijff(const RooCruijff &other, const char *name)
   : RooAbsPdf{other, name},
     _x{"x", this, other._x},
     _mu{"mean", this, other._mu},
     _sigmaL{"sigmaL", this, other._sigmaL},
     _sigmaR{"sigmaR", this, other._sigmaR},
     _alphaL{"alphaL", this, other._alphaL},
     _alphaR{"alphaR", this, other._alphaR},
     _beta{"beta", this, other._beta}
{
}

double RooCruijff::evaluate() const
{
   const double sigmaL = _sigmaL;
   const double sigmaR = _sigmaR;
   const double beta = _beta;
   const double arg = _x - _mu;
   const double arg2 = arg * arg;
   const double scale2 = 1.0 + 2 * beta * arg + beta * beta * arg2;
   return std::exp(arg < 0.0 ? -(arg2 * scale2) / (2 * sigmaL * sigmaL + _alphaL * arg2)
                             : -(arg2 * scale2) / (2 * sigmaR * sigmaR + _alphaR * arg2));
}

void RooCruijff::computeBatch(double *output, size_t nEvents, RooFit::Detail::DataMap const &dataMap) const
{
   RooBatchCompute::compute(dataMap.config(this), RooBatchCompute::Cruijff, output, nEvents,
                            {dataMap.at(_x), dataMap.at(_mu), dataMap.at(_sigmaL), dataMap.at(_sigmaR),
                             dataMap.at(_alphaL), dataMap.at(_alphaR), dataMap.at(_beta)});
}
