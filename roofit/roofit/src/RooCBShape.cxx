/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooCBShape
    \ingroup Roofit

PDF implementing the Crystal Ball line shape.
**/

#include "RooCBShape.h"

#include "RooRealVar.h"
#include "RooMath.h"
#include "RooBatchCompute.h"

#include <RooFit/Detail/MathFuncs.h>

#include "TMath.h"

#include <cmath>


////////////////////////////////////////////////////////////////////////////////

RooCBShape::RooCBShape(const char *name, const char *title,
             RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
             RooAbsReal& _alpha, RooAbsReal& _n) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alpha("alpha", "Alpha", this, _alpha),
  n("n", "Order", this, _n)
{
}

////////////////////////////////////////////////////////////////////////////////

RooCBShape::RooCBShape(const RooCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alpha("alpha", this, other.alpha),
  n("n", this, other.n)
{
}

////////////////////////////////////////////////////////////////////////////////

double RooCBShape::evaluate() const
{
   return RooFit::Detail::MathFuncs::cbShape(m, m0, sigma, alpha, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of Crystal ball Shape distribution.
void RooCBShape::doEval(RooFit::EvalContext &ctx) const
{
   RooBatchCompute::compute(ctx.config(this), RooBatchCompute::CBShape, ctx.output(),
                            {ctx.at(m), ctx.at(m0), ctx.at(sigma), ctx.at(alpha), ctx.at(n)});
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

double RooCBShape::analyticalIntegral(Int_t /*code*/, const char *rangeName) const
{
   using namespace RooFit::Detail::MathFuncs;
   return cbShapeIntegral(m.min(rangeName), m.max(rangeName), m0, sigma, alpha, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise that we know the maximum of self for given (m0,alpha,n,sigma)

Int_t RooCBShape::getMaxVal(const RooArgSet& vars) const
{
   RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
     return 1 ;
   }
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return an upper bound on the *unnormalized* value of this pdf, as required
/// by the RooAbsReal::getMaxVal() contract.
///
/// The Gaussian core is `exp(-0.5 * t^2)`, which is at most one, and the
/// power-law tail is attached continuously at `t = -|alpha|` where its value is
/// `exp(-0.5 * alpha^2) <= 1` and from where it falls off monotonically. The
/// maximum of the shape is therefore exactly one, reached at `m == m0`.

double RooCBShape::maxVal(Int_t code) const
{
  R__ASSERT(code==1) ;

  return 1.0 ;
}
