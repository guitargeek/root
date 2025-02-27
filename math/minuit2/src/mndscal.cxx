// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include <cstddef>

namespace ROOT {

namespace Minuit2 {

void Mndscal(unsigned int n, double da, double *dx)
{
   for (unsigned int i = 0; i < n; ++i) {
      dx[i] *= da;
   }
}

} // namespace Minuit2

} // namespace ROOT
