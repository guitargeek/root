// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_MnParabolaFactory
#define ROOT_MyMinuit2_MnParabolaFactory

namespace ROOT {

namespace MyMinuit2 {

class MnParabola;
class MnParabolaPoint;

class MnParabolaFactory {

public:
   MnParabolaFactory() {}

   ~MnParabolaFactory() {}

   MnParabola operator()(const MnParabolaPoint &, const MnParabolaPoint &, const MnParabolaPoint &) const;

   MnParabola operator()(const MnParabolaPoint &, double, const MnParabolaPoint &) const;

private:
};

} // namespace MyMinuit2

} // namespace ROOT

#endif // ROOT_MyMinuit2_MnParabolaFactory
