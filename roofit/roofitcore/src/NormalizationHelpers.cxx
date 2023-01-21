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

#include <unordered_map>

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
   if (dynamic_cast<RooSimultaneous const *>(&arg)) {
      RooFit::CompileContext ctx;
      std::unique_ptr<RooAbsArg> head = arg.compileForNormSet(normSet, ctx);
      auto *simPdf = static_cast<RooSimultaneous *>(head.get());

      RooArgList newServers;

      for (auto *cat : static_range_cast<RooAbsCategoryLValue *>(simPdf->flattenedCatList())) {

         for (auto const &catState : *cat) {
            std::string const &catName = catState.first;

            if (RooAbsPdf *pdf = simPdf->getPdf(catName.c_str())) {

               auto binnedInfo = RooHelpers::getBinnedL(*pdf);

               const std::string origname = pdf->GetName();
               pdf = binnedInfo.actualPdf ? binnedInfo.actualPdf : pdf;

               if (binnedInfo.isBinnedL) {
                  pdf->setAttribute("BinnedLikelihoodActive");
               }

               std::unique_ptr<RooArgSet> pdfNormSet(static_cast<RooArgSet *>(
                  std::unique_ptr<RooArgSet>(pdf->getVariables())->selectByAttrib("__obs__", true)));

               std::unique_ptr<RooAbsArg> pdfClone = compileForNormSetImpl(*pdf, *pdfNormSet);

               pdfClone->setAttribute(("ORIGNAME:" + origname).c_str());
               newServers.addOwned(std::move(pdfClone));
            }
         }
      }

      std::unique_ptr<RooAbsArg> indexCatClone = compileForNormSetImpl(simPdf->indexCat(), normSet);
      indexCatClone->setAttribute((std::string("ORIGNAME:") + indexCatClone->GetName()).c_str());
      newServers.addOwned(std::move(indexCatClone));

      simPdf->redirectServers(newServers, true, true);

      // This hack is necessary because the owned components can't contain two
      // args with the same name, as the container is a RooArgSet. We work
      // around this by letting the first owned components own the new owned
      // components.
      // const_cast<RooArgSet &>(*simPdf->ownedComponents())[0]->addOwnedComponents(std::move(newServers));
      newServers.releaseOwnership(); // INTENTIONAL LEAK FOR NOW!

      return head;
   } else {
      RooFit::CompileContext ctx;
      std::unique_ptr<RooAbsArg> head = arg.compileForNormSet(normSet, ctx);
      ctx.compileServers(*head, normSet);
      return head;
   }
}

} // namespace Detail

} // namespace RooFit
