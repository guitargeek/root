// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MatrixInverse
#define ROOT_MyMinuit2_MatrixInverse

#include "MyMinuit2/ABTypes.h"
#include "MyMinuit2/ABObj.h"

namespace ROOT {

namespace MyMinuit2 {

template <class mtype, class M, class T>
class MatrixInverse {

public:
   MatrixInverse(const M &obj) : fObject(obj) {}

   ~MatrixInverse() {}

   typedef mtype Type;

   const M &Obj() const { return fObject; }

private:
   M fObject;
};

template <class M, class T>
class MatrixInverse<vec, M, T> {

private:
   MatrixInverse(const M &obj) : fObject(obj) {}

public:
   ~MatrixInverse() {}

   typedef vec Type;

   const M &Obj() const { return fObject; }

private:
   M fObject;
};

template <class mt, class M, class T>
inline ABObj<mt, MatrixInverse<mt, ABObj<mt, M, T>, T>, T> Inverse(const ABObj<mt, M, T> &obj)
{
   return ABObj<mt, MatrixInverse<mt, ABObj<mt, M, T>, T>, T>(MatrixInverse<mt, ABObj<mt, M, T>, T>(obj));
}

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_MatrixInverse
