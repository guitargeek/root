// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnMatrix.h"

#include <numeric>

namespace ROOT {

namespace Minuit2 {

double similarity(const LAVector &avec, const LASymMatrix &mat)
{
   // calculate the similarity vector-matrix product: V^T M V
   // use matrix product and then dot function

   LAVector tmp = mat * avec;

   return std::inner_product(avec.Data(), avec.Data() + avec.size(), tmp.Data(), 0.0);
}

} // namespace Minuit2

} // namespace ROOT
