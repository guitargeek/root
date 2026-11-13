/// \cond ROOFIT_INTERNAL

/*
 * Project: RooFit
 *
 * Copyright (c) 2026, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_RooFitLegacy_RooExpensiveObjectCache_h
#define RooFit_RooFitLegacy_RooExpensiveObjectCache_h

#include <TObject.h>
#include <TString.h>

#include <map>

/// Only kept to read the RooWorkspaces written before ROOT 6.44, which have a
/// RooExpensiveObjectCache as data member. It has no functionality.
class RooExpensiveObjectCache : public TObject {
public:
   RooExpensiveObjectCache() = default;
   RooExpensiveObjectCache(const RooExpensiveObjectCache &) = delete;
   ~RooExpensiveObjectCache() override
   {
      for (auto &item : _map) {
         delete item.second;
      }
   }

   class ExpensiveObject {
   public:
      ExpensiveObject() = default;
      ExpensiveObject(const ExpensiveObject &) = delete;
      virtual ~ExpensiveObject() { delete _payload; }

   protected:
      Int_t _uid = 0;                           ///< Unique element ID ;
      TObject *_payload = nullptr;              ///< Payload
      std::map<TString, double> _realRefParams; ///< Names and values of real-valued reference parameters
      std::map<TString, Int_t> _catRefParams;   ///< Names and values of discrete-valued reference parameters
      TString _ownerName;                       ///< Name of RooAbsArg object that is associated to cache contents

      ClassDef(ExpensiveObject, 2);
   };

protected:
   Int_t _nextUID = 0;

   std::map<TString, ExpensiveObject *> _map;

   ClassDefOverride(RooExpensiveObjectCache, 2)
};

#endif

/// \endcond
