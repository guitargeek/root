// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Types                                                                 *
 *                                             *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

/*! \class TMVA::Types
\ingroup TMVA
Singleton class for Global types used by TMVA
*/

#include "TMVA/Types.h"

#include "TMVA/MsgLogger.h"

#include "RtypesCore.h"
#include "TString.h"
#include "TROOT.h"

#include <map>
#include <atomic>

namespace {

auto &Log()
{
   static TMVA::MsgLogger logger{"Types"};
   return logger;
}

auto &str2type()
{
   static std::map<TString, TMVA::Types::EMVA> map; // types-to-text map
   return map;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::Types::AddTypeMapping( Types::EMVA method, const TString& methodname )
{
   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   auto it = str2type().find( methodname );
   if (it != str2type().end()) {
      Log() << kFATAL
            << "Cannot add method " << methodname
            << " to the name->type map because it exists already" << Endl;
      return kFALSE;
   }

   str2type()[methodname] = method;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the method type (enum) for a given method (string)

TMVA::Types::EMVA TMVA::Types::GetMethodType( const TString& method )
{
   auto it = str2type().find( method );
   if (it == str2type().end()) {
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
      Log() << kFATAL << "Unknown method in map: " << method << Endl;
      return kVariable; // Inserted to get rid of GCC warning...
   }
   else return it->second;
}

////////////////////////////////////////////////////////////////////////////////

TString TMVA::Types::GetMethodName( TMVA::Types::EMVA method )
{
   auto it = str2type().begin();
   for (; it!=str2type().end(); ++it) if (it->second == method) return it->first;
   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
   Log() << kFATAL << "Unknown method index in map: " << method << Endl;
   return "";
}
