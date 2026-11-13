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
#ifndef ROO_OBJ_CACHE_MANAGER
#define ROO_OBJ_CACHE_MANAGER

#include "Rtypes.h"

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsCache.h"
#include "RooAbsCacheElement.h"
#include "RooCacheManager.h"


class RooObjCacheManager : public RooCacheManager<RooAbsCacheElement> {

public:

  RooObjCacheManager(RooAbsArg* owner=nullptr, Int_t maxSize=2, bool clearCacheOnServerRedirect=true, bool allowOptimize=false) ;
  RooObjCacheManager(const RooObjCacheManager& other, RooAbsArg* owner=nullptr) ;
  ~RooObjCacheManager() override ;

  bool redirectServersHook(const RooAbsCollection& /*newServerList*/, bool /*mustReplaceAll*/, bool /*nameChange*/, bool /*isRecursive*/) override ;
  void printCompactTreeHook(std::ostream&, const char *) override ;

  void insertObjectHook(RooAbsCacheElement&) override ;

  void setClearOnRedirect(bool flag) { _clearOnRedirect = flag ; }

protected:

  bool _clearOnRedirect ;
  bool _allowOptimize ; ///< No longer used, but kept because it is a persistent member (schema evolution)

  ClassDefOverride(RooObjCacheManager,3) ///< Cache manager for generic caches that contain RooAbsArg objects
} ;



#endif
