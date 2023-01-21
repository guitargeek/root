/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2022
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "NormalizationHelpers.h"

#include <RooAbsArg.h>
#include <RooArgList.h>
#include <RooHelpers.h>
#include <RooSimultaneous.h>

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

RooFit::CompileContext::~CompileContext() {}

void RooFit::CompileContext::compileServers(RooAbsArg &arg, RooArgSet const &normSet)
{
   auto serverClones = std::make_unique<RooArgList>();
   for (const auto server : arg.servers()) {
      auto serverClone = this->compile(*server, arg, normSet);
      if (serverClone) {
         serverClones->add(*serverClone);
      } else {
         this->compileServers(*server, normSet);
      }
   }
   arg.redirectServers(*serverClones, false, true);
}

RooAbsArg *RooFit::CompileContext::compile(RooAbsArg &arg, RooAbsArg &caller, RooArgSet const &normSet)
{
   if (auto existingServerClone = this->find(arg)) {
      return existingServerClone;
   }
   if (arg.isFundamental() && !normSet.find(arg)) {
      return nullptr;
   }
   if (arg.getAttribute("_COMPILED")) {
      return nullptr;
   }

   std::unique_ptr<RooAbsArg> newArg = arg.compileForNormSet(normSet, *this);
   this->compileServers(*newArg, normSet);
   const std::string attrib = std::string("ORIGNAME:") + arg.GetName();
   newArg->setAttribute(attrib.c_str());
   this->add(*newArg);
   RooAbsArg *out = newArg.get();
   caller.addOwnedComponents(std::move(newArg));
   return out;
}

namespace RooFit {

namespace Detail {

std::unique_ptr<RooAbsArg> compileForNormSetImpl(RooAbsArg const &arg, RooArgSet const &normSet)
{
   RooFit::CompileContext ctx;
   std::unique_ptr<RooAbsArg> head = arg.compileForNormSet(normSet, ctx);
   if (!dynamic_cast<RooSimultaneous const *>(&arg)) {
      ctx.compileServers(*head, normSet);
   }
   return head;
}

} // namespace Detail

} // namespace RooFit
