// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MnTiny
#define ROOT_MyMinuit2_MnTiny

namespace ROOT {

namespace MyMinuit2 {

class MnTiny {

public:
   MnTiny() : fOne(1.) {}

   ~MnTiny() {}

   double One() const;

   double operator()(double epsp1) const;

private:
   double fOne;
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_MnTiny
