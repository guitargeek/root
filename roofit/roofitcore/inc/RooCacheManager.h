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
#ifndef ROO_CACHE_MANAGER
#define ROO_CACHE_MANAGER

#include <RooAbsCache.h>

#include <Rtypes.h>

#include <vector>

class RooAbsCacheElement;
class RooArgSet;
class RooNormSetCache;

class RooCacheManager : public RooAbsCache {

public:

  RooCacheManager(int maxSize=2) ;
  RooCacheManager(RooAbsArg* owner, int maxSize=2) ;
  RooCacheManager(const RooCacheManager& other, RooAbsArg* owner=nullptr) ;
  ~RooCacheManager() override ;

  /// Getter function without integration set
  RooAbsCacheElement* getObj(const RooArgSet* nset, int* sterileIndex=nullptr, const TNamed* isetRangeName=nullptr) {
    return getObj(nset,nullptr,sterileIndex,isetRangeName) ;
  }

  /// Setter function without integration set
  int setObj(const RooArgSet* nset, RooAbsCacheElement* obj, const TNamed* isetRangeName=nullptr) {
    return setObj(nset,nullptr,obj,isetRangeName) ;
  }

  // Not worth to inline because it uses RooNameReg::ptr, which is expensive
  RooAbsCacheElement* getObj(const RooArgSet* nset, const RooArgSet* iset, int* sterileIdx, const char* isetRangeName);

  RooAbsCacheElement* getObj(const RooArgSet* nset, const RooArgSet* iset, int* sterileIndex=nullptr, const TNamed* isetRangeName=nullptr) ;
  int setObj(const RooArgSet* nset, const RooArgSet* iset, RooAbsCacheElement* obj, const TNamed* isetRangeName=nullptr) ;

  void reset() ;
  virtual void sterilize() ;

  /// Return index of slot used in last get or set operation
  int lastIndex() const {
    return _lastIndex ;
  }
  /// Return size of cache
  int cacheSize() const {
    return _size ;
  }

  /// Interface function to intercept server redirects
  bool redirectServersHook(const RooAbsCollection& /*newServerList*/, bool /*mustReplaceAll*/,
                 bool /*nameChange*/, bool /*isRecursive*/) override {
    return false ;
  }
  /// Interface function to intercept cache operation mode changes
  void operModeHook() override {
  }
  /// Interface function to cache add contents to output in tree printing mode
  void printCompactTreeHook(std::ostream&, const char *) override {
  }

  RooAbsCacheElement* getObjByIndex(int index) const ;
  RooArgSet selectFromSet1(RooArgSet const& argSet, int index) const ;
  RooArgSet selectFromSet2(RooArgSet const& argSet, int index) const ;

  /// Interface function to perform post-insert operations on cached object
  virtual void insertObjectHook(RooAbsCacheElement&) {
  }

  void wireCache() override;

protected:

  int _maxSize ;    ///<! Maximum size
  int _size ;       ///<! Actual use
  int _lastIndex ;  ///<! Last slot accessed

  std::vector<RooNormSetCache> _nsetCache ; ///<! Normalization/Integration set manager
  std::vector<RooAbsCacheElement*> _object ;                 ///<! Payload
  bool _wired ;               ///<! In wired mode, there is a single payload which is returned always
} ;

#endif
