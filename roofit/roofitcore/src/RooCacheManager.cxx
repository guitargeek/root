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
\file RooCacheManager.cxx
\class RooCacheManager
\ingroup Roofitcore

Template class RooCacheManager manages the storage of any type of data indexed on
the choice of normalization and optionally the set of integrated observables.
The purpose of this class is to faciliate storage of intermediate results
in operator p.d.f.s whose value and inner working are often highly dependent
on the user provided choice of normalization in getVal().

For efficiency reasons these normalization set pointer are
derefenced as little as possible. This class contains a lookup
table for RooArgSet pointer pairs -> normalization lists.  Distinct
pointer pairs that represent the same normalization/projection are
recognized and will all point to the same normalization list. Lists
for up to 'maxSize' different normalization/ projection
configurations can be cached.
**/

#include <RooCacheManager.h>

#include <RooAbsCacheElement.h>
#include <RooArgSet.h>
#include <RooHelpers.h>
#include <RooMsgService.h>
#include <RooNameReg.h>
#include <RooNormSetCache.h>

/// Constructor for simple caches without RooAbsArg payload. A cache
/// made with this constructor is not registered with its owner
/// and will not receive information on server redirects and
/// cache operation mode changes.
RooCacheManager::RooCacheManager(Int_t maxSize) : RooAbsCache(0)
{
   _maxSize = maxSize;
   _nsetCache.resize(_maxSize); // = new RooNormSetCache[maxSize] ;
   _object.resize(_maxSize, 0); // = new RooAbsCacheElement*[maxSize] ;
   _wired = false;
}

/// Constructor for simple caches with RooAbsArg derived payload. A cache
/// made with this constructor is registered with its owner
/// and will receive information on server redirects and
/// cache operation mode changes.
RooCacheManager::RooCacheManager(RooAbsArg *owner, Int_t maxSize) : RooAbsCache(owner)
{
   _maxSize = maxSize;
   _size = 0;

   _nsetCache.resize(_maxSize); // = new RooNormSetCache[maxSize] ;
   _object.resize(_maxSize, 0); // = new RooAbsCacheElement*[maxSize] ;
   _wired = false;
   _lastIndex = -1;

   Int_t i;
   for (i = 0; i < _maxSize; i++) {
      _object[i] = 0;
   }
}

/// Copy constructor.
RooCacheManager::RooCacheManager(const RooCacheManager &other, RooAbsArg *owner) : RooAbsCache(other, owner)
{
   _maxSize = other._maxSize;
   _size = other._size;

   _nsetCache.resize(_maxSize); // = new RooNormSetCache[_maxSize] ;
   _object.resize(_maxSize, 0); // = new RooAbsCacheElement*[_maxSize] ;
   _wired = false;
   _lastIndex = -1;

   // std::cout << "RooCacheManager:cctor(" << this << ") other = " << &other << " _size=" << _size << " _maxSize = " <<
   // _maxSize << std::endl ;

   Int_t i;
   for (i = 0; i < other._size; i++) {
      _nsetCache[i].initialize(other._nsetCache[i]);
      _object[i] = 0;
   }

   for (i = other._size; i < _maxSize; i++) {
      _object[i] = 0;
   }
}

/// Destructor
RooCacheManager::~RooCacheManager()
{
   for (int i = 0; i < _size; i++) {
      delete _object[i];
   }
}

/// Clear the cache
void RooCacheManager::reset()
{
   for (int i = 0; i < _maxSize; i++) {
      delete _object[i];
      _object[i] = 0;
      _nsetCache[i].clear();
   }
   _lastIndex = -1;
   _size = 0;
}

/// Clear the cache payload but retain slot mapping w.r.t to
/// normalization and integration sets.
void RooCacheManager::sterilize()
{
   Int_t i;
   for (i = 0; i < _maxSize; i++) {
      delete _object[i];
      _object[i] = 0;
   }
}

/// Insert payload object 'obj' in cache indexed on nset,iset and isetRangeName.
Int_t RooCacheManager::setObj(const RooArgSet *nset, const RooArgSet *iset, RooAbsCacheElement *obj,
                              const TNamed *isetRangeName)
{
   // Check if object is already registered
   Int_t sterileIdx(-1);
   if (getObj(nset, iset, &sterileIdx, isetRangeName)) {
      delete obj; // important! do not forget to cleanup memory
      return lastIndex();
   }

   if (sterileIdx >= 0) {
      // Found sterile slot that can should be recycled [ sterileIndex only set if isetRangeName matches ]

      if (sterileIdx >= _maxSize) {
         // cout << "RooCacheManager::setObj()/SI increasing object cache size from " << _maxSize << " to " <<
         // sterileIdx+4 << endl ;
         _maxSize = sterileIdx + 4;
         _object.resize(_maxSize, 0);
         _nsetCache.resize(_maxSize);
      }

      _object[sterileIdx] = obj;

      // Allow optional post-processing of object inserted in cache
      insertObjectHook(*obj);

      return lastIndex();
   }

   if (_size >= _maxSize - 1) {
      // cout << "RooCacheManager::setObj() increasing object cache size from " << _maxSize << " to " << _maxSize*2 <<
      // endl ;
      _maxSize *= 2;
      _object.resize(_maxSize, 0);
      _nsetCache.resize(_maxSize);
   }

   // cout << "RooCacheManager::setObj<T>(" << this << ") _size = " << _size << " _maxSize = " << _maxSize << endl ;
   _nsetCache[_size].autoCache(_owner, nset, iset, isetRangeName, true);
   if (_object[_size]) {
      delete _object[_size];
   }

   _object[_size] = obj;
   _size++;

   // Allow optional post-processing of object inserted in cache
   insertObjectHook(*obj);

   // Unwire cache in case it was wired
   _wired = false;

   return _size - 1;
}

RooAbsCacheElement *
RooCacheManager::getObj(const RooArgSet *nset, const RooArgSet *iset, Int_t *sterileIdx, const char *isetRangeName)
{
   // Not worth to inline because it uses RooNameReg::ptr, which is expensive
   if (_wired)
      return _object[0];
   return getObj(nset, iset, sterileIdx, RooNameReg::ptr(isetRangeName));
}

/// Retrieve payload object indexed on nset,uset amd isetRangeName
/// If sterileIdx is not null, it is set to the index of the sterile
/// slot in cacse such a slot is recycled.
RooAbsCacheElement *
RooCacheManager::getObj(const RooArgSet *nset, const RooArgSet *iset, Int_t *sterileIdx, const TNamed *isetRangeName)
{
   // Fast-track for wired mode
   if (_wired) {
      if (_object[0] == 0 && sterileIdx)
         *sterileIdx = 0;
      return _object[0];
   }

   Int_t i;
   for (i = 0; i < _size; i++) {
      if (_nsetCache[i].contains(nset, iset, isetRangeName) == true) {
         _lastIndex = i;
         if (_object[i] == 0 && sterileIdx)
            *sterileIdx = i;
         return _object[i];
      }
   }

   for (i = 0; i < _size; i++) {
      if (_nsetCache[i].autoCache(_owner, nset, iset, isetRangeName, false) == false) {
         _lastIndex = i;
         if (_object[i] == 0 && sterileIdx)
            *sterileIdx = i;
         return _object[i];
      }
   }

   return 0;
}

/// Retrieve payload object by slot index.
RooAbsCacheElement *RooCacheManager::getObjByIndex(Int_t index) const
{
   if (index < 0 || index >= _size) {
      oocoutE(_owner, ObjectHandling) << "RooCacheManager::getNormListByIndex: ERROR index (" << index
                                      << ") out of range [0," << _size - 1 << "]" << std::endl;
      return 0;
   }
   return _object[index];
}

/// Create RooArgSet contatining the objects that are both in the cached set 1
// with a given index and an input argSet.
RooArgSet RooCacheManager::selectFromSet1(RooArgSet const &argSet, int index) const
{
   return RooHelpers::selectFromArgSet(argSet, _nsetCache.at(index).nameSet1());
}

/// Create RooArgSet contatining the objects that are both in the cached set 2
// with a given index and an input argSet.
RooArgSet RooCacheManager::selectFromSet2(RooArgSet const &argSet, int index) const
{
   return RooHelpers::selectFromArgSet(argSet, _nsetCache.at(index).nameSet2());
}

void RooCacheManager::wireCache()
{
   if (_size == 0) {
      oocoutI(_owner, Optimization) << "RooCacheManager::wireCache(" << _owner->GetName() << ") no cached elements!"
                                    << std::endl;
   } else if (_size == 1) {
      oocoutI(_owner, Optimization) << "RooCacheManager::wireCache(" << _owner->GetName() << ") now wiring cache"
                                    << std::endl;
      _wired = true;
   } else if (_size > 1) {
      oocoutI(_owner, Optimization) << "RooCacheManager::wireCache(" << _owner->GetName()
                                    << ") cache cannot be wired because it contains more than one element" << std::endl;
   }
}
