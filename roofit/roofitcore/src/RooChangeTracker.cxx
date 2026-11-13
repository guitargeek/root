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
\file RooChangeTracker.cxx
\class RooChangeTracker
\ingroup Roofitcore

Meta object that tracks value changes in a given set of RooAbsArgs by
comparing their values with reference values.
**/

#include "RooChangeTracker.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooArgSet.h"


////////////////////////////////////////////////////////////////////////////////
/// Constructor. The set trackSet contains the arguments to be tracked for
/// changes.

RooChangeTracker::RooChangeTracker(const char* name, const char* title, const RooArgSet& trackSet) :
  RooAbsReal(name, title),
  _realSet("realSet","Set of real-valued components to be tracked",this),
  _catSet("catSet","Set of discrete-valued components to be tracked",this)
{
  for (const auto arg : trackSet) {
    if (dynamic_cast<const RooAbsReal*>(arg)) {
      _realSet.add(*arg) ;
    }
    if (dynamic_cast<const RooAbsCategory*>(arg)) {
      _catSet.add(*arg) ;
    }
  }

  _realRef.reserve(_realSet.size());
  for (auto *real : static_range_cast<RooAbsReal *>(_realSet)) {
    _realRef.push_back(real->getVal());
  }

  _catRef.reserve(_catSet.size());
  for (auto *cat : static_range_cast<RooAbsCategory *>(_catSet)) {
    _catRef.push_back(cat->getCurrentIndex());
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooChangeTracker::RooChangeTracker(const RooChangeTracker& other, const char* name) :
  RooAbsReal(other, name),
  _realSet("realSet",this,other._realSet),
  _catSet("catSet",this,other._catSet),
  _realRef(other._realRef),
  _catRef(other._catRef)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Returns true if a tracked value differs from its reference value. If
/// clearState is true, the current values become the new reference values.
/// The first call with clearState=true always returns true.

bool RooChangeTracker::hasChanged(bool clearState)
{
  bool valuesChanged = false;

  for (std::size_t i = 0; i < _realSet.size(); ++i) {
    const double val = static_cast<const RooAbsReal &>(_realSet[i]).getVal();
    if (val != _realRef[i]) {
      valuesChanged = true;
      if (clearState) _realRef[i] = val;
    }
  }

  for (std::size_t i = 0; i < _catSet.size(); ++i) {
    const Int_t idx = static_cast<const RooAbsCategory &>(_catSet[i]).getCurrentIndex();
    if (idx != _catRef[i]) {
      valuesChanged = true;
      if (clearState) _catRef[i] = idx;
    }
  }

  if (clearState && !_init) {
    valuesChanged = true;
    _init = true;
  }

  return valuesChanged;
}

////////////////////////////////////////////////////////////////////////////////

RooArgSet RooChangeTracker::parameters() const
{
  RooArgSet ret ;
  ret.add(_realSet) ;
  ret.add(_catSet) ;
  return ret ;
}
