// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/BasicMinimumError.h"

#include "MyMinuit2/MnMatrix.h"

#include "MyMinuit2/MnPrint.h"

namespace ROOT {

namespace MyMinuit2 {

MnAlgebraicSymMatrix BasicMinimumError::Hessian() const
{
   // calculate Heassian: inverse of error matrix
   MnAlgebraicSymMatrix tmp(fMatrix);
   int ifail = Invert(tmp);
   if (ifail != 0) {
      MnPrint print("BasicMinimumError::Hessian");
      print.Warn("Inversion fails; return diagonal matrix");
      MnAlgebraicSymMatrix tmp2(fMatrix.Nrow());
      for (unsigned int i = 0; i < fMatrix.Nrow(); i++) {
         tmp2(i, i) = 1. / fMatrix(i, i);
      }
      return tmp2;
   }
   return tmp;
}

} // namespace MyMinuit2

} // namespace ROOT
