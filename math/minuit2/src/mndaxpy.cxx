// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

namespace ROOT {

namespace Minuit2 {

void Mndaxpy(unsigned int n, double da, const double *dx, double *dy)
{
   for (unsigned int i = 0; i < n; ++i) {
      dy[i] += da * dx[i];
   }
}

} // namespace Minuit2

} // namespace ROOT
