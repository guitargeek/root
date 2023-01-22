/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2023
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <NormalizationHelpers.h>

#include <RooAbsArg.h>
#include <RooArgList.h>
#include <RooArgSet.h>

RooFit::CompileContext::CompileContext(RooArgSet const &topLevelNormSet) : _topLevelNormSet{topLevelNormSet} {}

RooFit::CompileContext::~CompileContext() {}

void RooFit::CompileContext::add(RooAbsArg &arg)
{
   _clonedArgsSet.emplace(arg.namePtr(), &arg);
}

RooAbsArg *RooFit::CompileContext::find(RooAbsArg &arg) const
{
   auto existingServerClone = _clonedArgsSet.find(arg.namePtr());
   if (existingServerClone != _clonedArgsSet.end()) {
      return existingServerClone->second;
   }
   return nullptr;
}

void RooFit::CompileContext::compileServers(RooAbsArg &arg, RooArgSet const &normSet)
{
   RooArgList serverClones;
   for (const auto server : arg.servers()) {
      if (auto serverClone = this->compile(*server, arg, normSet)) {
         serverClones.add(*serverClone);
      }
   }
   arg.redirectServers(serverClones, false, true);
}

RooAbsArg *RooFit::CompileContext::compile(RooAbsArg &arg, RooAbsArg &owner, RooArgSet const &normSet)
{
   if (auto existingServerClone = this->find(arg)) {
      return existingServerClone;
   }
   if (arg.isFundamental() && !_topLevelNormSet.find(arg)) {
      return nullptr;
   }
   if (arg.getAttribute("_COMPILED")) {
      return nullptr;
   }

   std::unique_ptr<RooAbsArg> newArg = arg.compileForNormSet(normSet, *this);
   newArg->setAttribute("_COMPILED");
   const std::string attrib = std::string("ORIGNAME:") + arg.GetName();
   newArg->setAttribute(attrib.c_str());
   this->add(*newArg);
   RooAbsArg *out = newArg.get();
   owner.addOwnedComponents(std::move(newArg));
   return out;
}
