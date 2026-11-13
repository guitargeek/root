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

/**
\file RooObjCacheManager.cxx
\class RooObjCacheManager
\ingroup Roofitcore

Implementation of a RooCacheManager<RooAbsCacheElement>
that specializes in the storage of cache elements that contain RooAbsArg objects.
Caches with RooAbsArg derived payload require special care as server redirects
may need to be forwarded to such cache payload. This cache manager takes care of this
by forwarding these calls to the RooAbsCacheElement interface functions, which
have a sensible default implementation.
**/

#include "RooObjCacheManager.h"
#include "RooMsgService.h"

#include <ostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
/// Constructor of object cache manager for given owner. If clearCacheOnServerRedirect is true
/// all cache elements will be cleared when a server redirect is intercepted by the cache manager.
/// This is the default strategy and should only be overridden when you really understand
/// what you're doing as properly implementing server redirect in cache elements can get very
/// complicated, especially if there are (cyclical) reference back to the owning object

RooObjCacheManager::RooObjCacheManager(RooAbsArg* owner, Int_t maxSize, bool clearCacheOnServerRedirect, bool allowOptimize) :
  RooCacheManager<RooAbsCacheElement>(owner,maxSize),
  _clearOnRedirect(clearCacheOnServerRedirect),
  _allowOptimize(allowOptimize)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooObjCacheManager::RooObjCacheManager(const RooObjCacheManager& other, RooAbsArg* owner) :
  RooCacheManager<RooAbsCacheElement>(other,owner),
  _clearOnRedirect(other._clearOnRedirect),
  _allowOptimize(other._allowOptimize)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooObjCacheManager::~RooObjCacheManager() = default;


////////////////////////////////////////////////////////////////////////////////
/// Intercept server redirect calls. If clearOnRedirect was set, sterilize
/// the cache (i.e. keep the structure but delete all contents). If not
/// forward serverRedirect to cache elements

bool RooObjCacheManager::redirectServersHook(const RooAbsCollection& newServerList, bool mustReplaceAll, bool nameChange, bool isRecursive)
{
  if (_clearOnRedirect) {

    sterilize() ;

  } else {

    for (Int_t i=0 ; i<cacheSize() ; i++) {
      if (_object[i]) {
   _object[i]->redirectServersHook(newServerList,mustReplaceAll,nameChange,isRecursive) ;
      }
    }

  }

  return false ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set owner link on all object inserted into cache.

void RooObjCacheManager::insertObjectHook(RooAbsCacheElement& obj)
{
  obj.setOwner(_owner) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Add details on cache contents when printing in tree mode

void RooObjCacheManager::printCompactTreeHook(std::ostream& os, const char *indent)
{
  for (Int_t i=0 ; i<cacheSize() ; i++) {
    if (_object[i]) {
      _object[i]->printCompactTreeHook(os,indent,i,cacheSize()-1) ;
    }
  }
}


