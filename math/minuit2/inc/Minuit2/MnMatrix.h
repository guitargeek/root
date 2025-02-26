// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnMatrix
#define ROOT_Minuit2_MnMatrix

// add MnConfig file to define before everything compiler
// dependent macros

#include "Minuit2/MnConfig.h"

// define typedf's in MnMatrixfwd
#include "Minuit2/MnMatrixfwd.h"

#include "Minuit2/ABObj.h"
#include "Minuit2/LASymMatrix.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/MatrixInverse.h"

namespace ROOT {

namespace Minuit2 {

// Matrix-vector product
inline ABObj<vec, ABProd<ABObj<sym, LASymMatrix>, ABObj<vec, LAVector>>>
operator*(const ABObj<sym, LASymMatrix> &a, const ABObj<vec, LAVector> &b)
{
   return {ABProd<ABObj<sym, LASymMatrix>, ABObj<vec, LAVector>>(a, b)};
}

///    LAPACK Algebra functions
///    specialize the Invert function for LASymMatrix

inline ABObj<sym, MatrixInverse<sym, ABObj<sym, LASymMatrix>, double>>
Inverse(const ABObj<sym, LASymMatrix> &obj)
{
   return {MatrixInverse<sym, ABObj<sym, LASymMatrix>, double>{obj}};
}

int Invert(LASymMatrix &);

int Invert_undef_sym(LASymMatrix &);

///    LAPACK Algebra function
///    specialize the Outer_product function for LAVector;

inline ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector>, double>> Outer_product(const ABObj<vec, LAVector> &obj)
{
   return {VectorOuterProduct<ABObj<vec, LAVector>, double>{obj}};
}

void Outer_prod(LASymMatrix &, const LAVector &, double f = 1.);

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnMatrix
