/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

#include "Riostream.h"
#include "RooAbsProxy.h"
#include "RooArgSet.h"

/**
\file RooAbsProxy.cxx
\class RooAbsProxy
\ingroup Roofitcore

Abstract interface for proxy classes.
Proxy classes hold pointers to other Roofit objects
and process serverRedirect changes so that the proxied
pointers are updated accordingly on a clone or copy of
of the owning class
**/



RooAbsProxy::RooAbsProxy() = default;

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsProxy::RooAbsProxy(const char* /*name*/, const RooAbsProxy& other) :
  _nset(other._nset)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

void RooAbsProxy::changeNormSet(const RooArgSet* newNormSet)
{
  _nset = const_cast<RooArgSet*>(newNormSet) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print proxy name

void RooAbsProxy::print(std::ostream& os, bool /*addContents*/) const
{
  os << name() << std::endl;
}
