// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_DavidonErrorUpdator
#define ROOT_MyMinuit2_DavidonErrorUpdator

#include "MyMinuit2/MinimumErrorUpdator.h"

namespace ROOT {

   namespace MyMinuit2 {


/**
   Update of the covariance matrix for the Variable Metric minimizer (MIGRAD)
 */
class DavidonErrorUpdator : public MinimumErrorUpdator {

public:

  DavidonErrorUpdator() {}

  virtual ~DavidonErrorUpdator() {}

  virtual MinimumError Update(const MinimumState&, const MinimumParameters&,
                              const FunctionGradient&) const;

private:

};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_DavidonErrorUpdator
