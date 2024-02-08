// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_VectorOuterProduct
#define ROOT_MyMinuit2_VectorOuterProduct

#include "MyMinuit2/ABTypes.h"
#include "MyMinuit2/ABObj.h"

namespace ROOT {

namespace MyMinuit2 {

template <class M, class T>
class VectorOuterProduct {

public:
   VectorOuterProduct(const M &obj) : fObject(obj) {}

   ~VectorOuterProduct() {}

   typedef sym Type;

   const M &Obj() const { return fObject; }

private:
   M fObject;
};

template <class M, class T>
inline ABObj<sym, VectorOuterProduct<ABObj<vec, M, T>, T>, T> Outer_product(const ABObj<vec, M, T> &obj)
{
   return ABObj<sym, VectorOuterProduct<ABObj<vec, M, T>, T>, T>(VectorOuterProduct<ABObj<vec, M, T>, T>(obj));
}

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_VectorOuterProduct
