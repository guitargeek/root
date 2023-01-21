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

namespace {

void addServerClonesToList(RooAbsArg &arg, std::unordered_map<TNamed const *, RooAbsArg *> &clonedArgsSet,
                           RooArgSet const &normSet)
{
   RooArgList serverClones;
   RooArgList owned;
   for (const auto server : arg.servers()) {
      auto existingServerClone = clonedArgsSet.find(server->namePtr());
      if (existingServerClone != clonedArgsSet.end()) {
         serverClones.add(*existingServerClone->second);
      } else if (!server->isFundamental() || normSet.find(*server)) {

         if (!server->getAttribute("_COMPILED")) {

            std::unique_ptr<RooAbsArg> serverClone = server->compileForNormSet(normSet, normSet);
            addServerClonesToList(*serverClone, clonedArgsSet, normSet);
            const std::string attrib = std::string("ORIGNAME:") + server->GetName();
            serverClone->setAttribute(attrib.c_str());
            clonedArgsSet.emplace(serverClone->namePtr(), serverClone.get());
            serverClones.add(*serverClone);
            owned.addOwned(std::move(serverClone));

         } else {
            addServerClonesToList(*server, clonedArgsSet, normSet);
         }
      }
   }
   arg.redirectServers(serverClones, false, true);
   arg.addOwnedComponents(std::move(owned));
}

} // namespace

namespace RooFit {

namespace Detail {

std::unique_ptr<RooAbsArg> compileForNormSetImpl(RooAbsArg const &arg, RooArgSet const &normSet)
{
   std::unique_ptr<RooAbsArg> head = arg.compileForNormSet(normSet, normSet);

   if (auto *simPdf = dynamic_cast<RooSimultaneous *>(head.get())) {
      for (auto *cat : static_range_cast<RooAbsCategoryLValue *>(simPdf->flattenedCatList())) {

         for (auto const &catState : *cat) {
            std::string const &catName = catState.first;

            RooArgList newServers;
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
            simPdf->redirectServers(newServers, false, true);

            simPdf->addOwnedComponents(std::move(newServers));
         }
      }
   } else {
      std::unordered_map<TNamed const *, RooAbsArg *> clonedArgsSet;
      addServerClonesToList(*head, clonedArgsSet, normSet);
   }

   return head;
}

} // namespace Detail

} // namespace RooFit
