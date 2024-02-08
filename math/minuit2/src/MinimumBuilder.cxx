// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "MyMinuit2/MinimumBuilder.h"
#include "MyMinuit2/MnPrint.h"

namespace ROOT {

namespace MyMinuit2 {

MinimumBuilder::MinimumBuilder() : fPrintLevel(MnPrint::GlobalLevel()), fStorageLevel(1), fTracer(nullptr) {}

} // namespace MyMinuit2

} // namespace ROOT
