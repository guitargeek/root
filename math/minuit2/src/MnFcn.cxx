// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/MnFcn.h"
#include "MyMinuit2/FCNBase.h"
#include "MyMinuit2/MnVectorTransform.h"

namespace ROOT {

   namespace MyMinuit2 {


MnFcn::~MnFcn() {
   //   std::cout<<"Total number of calls to FCN: "<<fNumCall<<std::endl;
}

double MnFcn::operator()(const MnAlgebraicVector& v) const {
   // evaluate FCN converting from from MnAlgebraicVector to std::vector
   fNumCall++;
   return fFCN(MnVectorTransform()(v));
}

// double MnFcn::operator()(const std::vector<double>& par) const {
//     return fFCN(par);
// }

double MnFcn::ErrorDef() const {return fFCN.Up();}

double MnFcn::Up() const {return fFCN.Up();}

   }  // namespace MyMinuit2

}  // namespace ROOT
