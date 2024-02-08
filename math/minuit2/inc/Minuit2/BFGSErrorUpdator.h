// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_MyMinuit2_BFGSErrorUpdator
#define ROOT_MyMinuit2_BFGSErrorUpdator

#include "MyMinuit2/MinimumErrorUpdator.h"

namespace ROOT {

   namespace MyMinuit2 {


/**
   Update of the covariance matrix for the Variable Metric minimizer (MIGRAD)
 */
class BFGSErrorUpdator : public MinimumErrorUpdator {

public:

  BFGSErrorUpdator() {}

  virtual ~BFGSErrorUpdator() {}

  virtual MinimumError Update(const MinimumState&, const MinimumParameters&,
                              const FunctionGradient&) const;

private:

};

  }  // namespace MyMinuit2

}  // namespace ROOT

#endif  // ROOT_MyMinuit2_BFGSErrorUpdator
