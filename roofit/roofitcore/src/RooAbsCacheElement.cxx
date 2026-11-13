/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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

/**
\file RooAbsCacheElement.cxx
\class RooAbsCacheElement
\ingroup Roofitcore

Abstract base class for objects to be stored
in RooAbsCache cache manager objects. Each storage element has an
interface to pass on calls for server redirection.
**/

#include "RooAbsCacheElement.h"
#include "RooAbsArg.h"
#include "RooArgList.h"



////////////////////////////////////////////////////////////////////////////////
/// Interface for server redirect calls

bool RooAbsCacheElement::redirectServersHook(const RooAbsCollection& /*newServerList*/, bool /*mustReplaceAll*/,
                      bool /*nameChange*/, bool /*isRecursive*/)
{
  return false ;
}



////////////////////////////////////////////////////////////////////////////////
/// Hook function to print cache guts in tree printing mode of RooAbsArgs

void RooAbsCacheElement::printCompactTreeHook(std::ostream&, const char *, Int_t , Int_t )
{
}


