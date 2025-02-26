// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ABObj
#define ROOT_Minuit2_ABObj

#include "Minuit2/ABTypes.h"

namespace ROOT {

namespace Minuit2 {

template <class Type, class M, class T>
class ABObj {

public:
   ABObj(const M &obj, T factor = 1.) : fObject(obj), fFactor(factor) {}

   // Assignment is historically deleted (that was done probably to avoid
   // mistakes from accidental re-assignment). Also define destructor and copy
   // constructor in this case, according to the rule of five.
   ~ABObj() = default;
   ABObj(const ABObj &) = default;
   ABObj(ABObj &&) = default;
   ABObj &operator=(const ABObj &) = delete;
   ABObj &operator=(ABObj &&) = delete;

   const M &Obj() const { return fObject; }

   T f() const { return fFactor; }

private:
   M fObject;
   T fFactor;
};

// templated scaling operator *
template <class mt, class M, class T>
ABObj<mt, M, T> operator*(T f, const M &obj)
{
   return {obj, f};
}

// templated operator /
template <class mt, class M, class T>
ABObj<mt, M, T> operator/(const M &obj, T f)
{
   return {obj, T(1.) / f};
}

// templated unary operator -
template <class mt, class M, class T>
ABObj<mt, M, T> operator-(const M &obj)
{
   return {obj, -1.};
}

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ABObj
