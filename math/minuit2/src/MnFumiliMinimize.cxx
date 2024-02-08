// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/MnFumiliMinimize.h"
#include "MyMinuit2/FunctionMinimum.h"
#include "MyMinuit2/FumiliMinimizer.h"

namespace ROOT {

   namespace MyMinuit2 {





FunctionMinimum MnFumiliMinimize::operator()(unsigned int maxfcn, double toler) {
   // minimize using Fumili
   // need to reimplement otherwise base class method is done

   assert(fState.IsValid());
   unsigned int npar = VariableParameters();
   //   assert(npar > 0);
   if(maxfcn == 0) maxfcn = 200 + 100*npar + 5*npar*npar;
   FunctionMinimum min = Minimizer().Minimize( Fcnbase(), fState, fStrategy, maxfcn, toler);
   fNumCall += min.NFcn();
   fState = min.UserState();
   return min;
}

   }  // namespace MyMinuit2

}  // namespace ROOT
