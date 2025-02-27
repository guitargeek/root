// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/LAVector.h"

#include <numeric>

namespace ROOT {

namespace Minuit2 {

double inner_product(const LAVector &v1, const LAVector &v2)
{
   // calculate inner (dot) product of two vectors
   return std::inner_product(v1.Data(), v1.Data() + v1.size(), v2.Data(), 0.0);
}

} // namespace Minuit2

} // namespace ROOT
