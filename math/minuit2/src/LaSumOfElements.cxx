// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"

#include <numeric>

namespace ROOT {

namespace Minuit2 {

double sum_of_elements(const LAVector &v)
{
   // calculate the absolute sum of the vector elements
   return std::accumulate(v.Data(), v.Data() + v.size(), 0.0);
}

double sum_of_elements(const LASymMatrix &m)
{
   // calculate the absolute sum of all the matrix elements
   return std::accumulate(m.Data(), m.Data() + m.size(), 0.0);
}

} // namespace Minuit2

} // namespace ROOT
