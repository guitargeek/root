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

/// \class RooMappedCategory
/// RooMappedCategory provides a category-to-category mapping defined
/// by pattern matching on their state labels.
///
/// The mapping function consists of a series of wild card regular expressions.
/// Each expression is matched to the input categories' state labels, and an associated
/// output state label.

#include "RooMappedCategory.h"

#include "RooStreamParser.h"
#include "RooMsgService.h"

#include "TBuffer.h"
#include "TString.h"
#include "TRegexp.h"

#include <ostream>

RooMappedCategory::RooMappedCategory() : _defCat(0) { }

RooMappedCategory::RooMappedCategory(const char *name, const char *title, RooAbsCategory& inputCat, const char* defOut, Int_t defOutIdx) :
  RooAbsCategory(name, title), _inputCat("input","Input category",this,inputCat)
{
  // Constructor with input category and name of default output state, which is assigned
  // to all input category states that do not follow any mapping rule.
  if (defOutIdx==NoCatIdx) {
    _defCat = defineState(defOut).second;
  } else {
    _defCat = defineState(defOut,defOutIdx).second;
  }
}

RooMappedCategory::RooMappedCategory(const RooMappedCategory &other, const char *name)
   : RooAbsCategory(other, name),
     _inputCat("input", this, other._inputCat),
     _mapArray(other._mapArray)
{
   setDefCat(lookupIndex(other.lookupName(other._defCat)));
}

RooMappedCategory::~RooMappedCategory() = default;

bool RooMappedCategory::map(const char* inKeyRegExp, const char* outKey, Int_t outIdx)
{
  // Add mapping rule: any input category state label matching the 'inKeyRegExp'
  // wildcard expression will be mapped to an output state with name 'outKey'
  //
  // Rules are evaluated in the order they were added. In case an input state
  // matches more than one rule, the first rules output state will be assigned

  if (!inKeyRegExp || !outKey) return true ;

  // Check if pattern is already registered
  if (_mapArray.find(inKeyRegExp)!=_mapArray.end()) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName() << "): ERROR expression "
                          << inKeyRegExp << " already mapped" << std::endl ;
    return true ;
  }

  // Check if output type exists, if not register
  value_type catIdx = lookupIndex(outKey);
  if (catIdx == invalidCategory().second) {
    if (outIdx==NoCatIdx) {
      catIdx = defineState(outKey).second;
    } else {
      catIdx = defineState(outKey,outIdx).second;
    }
  }

  if (catIdx == invalidCategory().second) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName()
                          << "): ERROR, unable to define category for output type " << outKey << std::endl ;
    return true;
  }

  // Create new map entry ;
  Entry e(inKeyRegExp, catIdx);
  if (!e.ok()) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName()
                          << "): ERROR, expression " << inKeyRegExp << " didn't compile" << std::endl ;
    return true ;
  }

  _mapArray[inKeyRegExp] = e ;
  return false ;
}



/// Return the output state of the first mapping rule that matches the current
/// label of the input category, or the default state if no rule matches.
RooAbsCategory::value_type RooMappedCategory::evaluate() const
{
   const char *inKey = _inputCat->getCurrentLabel();
   for (const auto &strAndEntry : _mapArray) {
      if (strAndEntry.second.match(inKey)) {
         return strAndEntry.second.outCat();
      }
   }
   return _defCat;
}

void RooMappedCategory::printMultiline(std::ostream& os, Int_t content, bool verbose, TString indent) const
{
  // Print info about this mapped category to the specified stream. In addition to the info
  // from RooAbsCategory::printStream() we add:
  //
  //  Standard : input category
  //     Shape : default value
  //   Verbose : list of mapping rules

  RooAbsCategory::printMultiline(os,content,verbose,indent);

  if (verbose) {
    os << indent << "--- RooMappedCategory ---" << std::endl
       << indent << "  Maps from " ;
    _inputCat.arg().printStream(os,0,kStandard);

    os << indent << "  Default value is " << lookupName(_defCat) << " = " << _defCat << '\n';

    os << indent << "  Mapping rules:" << std::endl;
    for (const auto& strAndEntry : _mapArray) {
      os << indent << "  " << strAndEntry.first << " -> " << strAndEntry.second.outCat() << std::endl ;
    }
  }
}


bool RooMappedCategory::readFromStream(std::istream& is, bool compact, bool /*verbose*/)
{
  // Read object contents from given stream
   if (compact) {
     coutE(InputArguments) << "RooMappedCategory::readFromSteam(" << GetName() << "): can't read in compact mode" << std::endl ;
     return true ;
   } else {

     //Clear existing definitions, but preserve default output
     TString defCatName(lookupName(_defCat));
     _mapArray.clear() ;
     clearTypes() ;
     _defCat = defineState(defCatName.Data()).second;

     TString token;
     TString errorPrefix("RooMappedCategory::readFromStream(");
     errorPrefix.Append(GetName()) ;
     errorPrefix.Append(")") ;
     RooStreamParser parser(is,errorPrefix) ;
     parser.setPunctuation(":,") ;

     TString destKey;
     TString srcKey;
     bool readToken(true) ;

    // Loop over definition sequences
     while(true) {
       if (readToken) token=parser.readToken() ;
       if (token.IsNull()) break ;
       readToken=true ;

       destKey = token ;
       if (parser.expectToken(":",true)) return true ;

       // Loop over list of sources for this destination
       while(true) {
         srcKey = parser.readToken() ;
         token = parser.readToken() ;

         // Map a value
         if (map(srcKey,destKey)) return true ;

         // Unless next token is ',' current token
         // is destination part of next sequence
         if (token.CompareTo(",")) {
           readToken = false ;
           break ;
         }
       }
     }
     return false ;
   }
   //return false ; // statement unreachable (OSF)
}


////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooMappedCategory to more intuitively reflect the contents of the
/// product operator construction

void RooMappedCategory::printMetaArgs(std::ostream& os) const
{
  // Scan array of regexps
  RooAbsCategory::value_type prevOutCat = invalidCategory().second;
  bool first(true) ;
  os << "map=(" ;
  for (const auto& iter : _mapArray) {
    if (iter.second.outCat() != prevOutCat) {
      if (!first) { os << " " ; }
      first=false ;

      os << iter.second.outCat() << ":" << iter.first ;
      prevOutCat = iter.second.outCat();
    } else {
      os << "," << iter.first ;
    }
  }

  if (!first) { os << " " ; }
  os << lookupName(_defCat) << ":*" ;

  os << ") " ;
}




void RooMappedCategory::writeToStream(std::ostream& os, bool compact) const
{
  // Write object contents to given stream
  if (compact) {
    // Write value only
    os << getCurrentLabel() ;
  } else {
    // Write mapping expression

    // Scan array of regexps
    RooAbsCategory::value_type prevOutCat = invalidCategory().second;
    bool first(true) ;
    for (const auto& iter : _mapArray) {
      if (iter.second.outCat() != prevOutCat) {
        if (!first) { os << " " ; }
        first=false ;

        os << iter.second.outCat() << "<-" << iter.first ;
        prevOutCat = iter.second.outCat();
      } else {
        os << "," << iter.first ;
      }
    }

    if (!first) { os << " " ; }
    os << lookupName(_defCat) << ":*" ;
  }
}


RooMappedCategory::Entry::Entry(const char *exp, RooAbsCategory::value_type cat) : _expr(exp), _catIdx(cat) {}
RooMappedCategory::Entry::Entry(const Entry &other) : _expr(other._expr), _catIdx(other._catIdx) {}
RooMappedCategory::Entry::~Entry() { delete _regexp; }
bool RooMappedCategory::Entry::ok() { return (const_cast<TRegexp*>(regexp())->Status()==TRegexp::kOK) ; }

////////////////////////////////////////////////////////////////////////////////

RooMappedCategory::Entry& RooMappedCategory::Entry::operator=(const RooMappedCategory::Entry& other)
{
  if (&other==this) return *this ;

  _expr = other._expr ;
  _catIdx = other._catIdx;

  if (_regexp) {
    delete _regexp ;
    _regexp = nullptr;
  }

  return *this;
}

bool RooMappedCategory::Entry::match(const char* testPattern) const {
  return (TString(testPattern).Index(*regexp())>=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Mangle name : escape regexp character '+'

TString RooMappedCategory::Entry::mangle(const char* exp) const
{
  TString t ;
  const char *c = exp ;
  while(*c) {
    if (*c=='+') t.Append('\\') ;
    t.Append(*c) ;
    c++ ;
  }
  return t ;
}

const TRegexp* RooMappedCategory::Entry::regexp() const {
  if (!_regexp) {
    _regexp = new TRegexp(mangle(_expr), true);
  }

  return _regexp;
}
