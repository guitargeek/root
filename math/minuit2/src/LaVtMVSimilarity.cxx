// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/LASymMatrix.h"
#include "MyMinuit2/LAVector.h"
#include "MyMinuit2/LaProd.h"

namespace ROOT {

   namespace MyMinuit2 {


double mnddot(unsigned int, const double*, int, const double*, int);

double similarity(const LAVector& avec, const LASymMatrix& mat) {
   // calculate the similarity vector-matrix product: V^T M V
   // use matrix product and then dot function (using mnddot)

   LAVector tmp = mat*avec;

   double value = mnddot(avec.size(), avec.Data(), 1, tmp.Data(), 1);
   return value;
}

   }  // namespace MyMinuit2

}  // namespace ROOT
