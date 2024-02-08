// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_FumiliGradientCalculator
#define ROOT_MyMinuit2_FumiliGradientCalculator

#include "MyMinuit2/GradientCalculator.h"
#include "MyMinuit2/MnMatrix.h"

namespace ROOT {

   namespace MyMinuit2 {



class FumiliFCNBase;
class MnUserTransformation;

class FumiliGradientCalculator : public GradientCalculator {

public:

  FumiliGradientCalculator(const FumiliFCNBase& fcn, const MnUserTransformation& state, int n) :
    fFcn(fcn),
    fTransformation(state),
    fHessian(MnAlgebraicSymMatrix(n) )
  {}

  ~FumiliGradientCalculator() {}

  FunctionGradient operator()(const MinimumParameters&) const;

  FunctionGradient operator()(const MinimumParameters&,
                              const FunctionGradient&) const;

  const MnUserTransformation& Trafo() const {return fTransformation;}

  const MnAlgebraicSymMatrix & Hessian() const { return fHessian; }

private:

  const FumiliFCNBase& fFcn;
  const MnUserTransformation& fTransformation;
  mutable MnAlgebraicSymMatrix fHessian;

};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_FumiliGradientCalculator
