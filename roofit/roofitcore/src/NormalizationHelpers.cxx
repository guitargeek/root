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

#include <RooAbsCachedPdf.h>
#include <RooAbsPdf.h>
#include <RooAddition.h>
#include <RooProdPdf.h>

#include "RooNormalizedPdf.h"

#include "RooBatchCompute.h"

/// A RooProdPdf with a fixed normalization set can be replaced by this class.
/// Its purpose is to provide the right client-server interface for the
/// evaluation of RooProdPdf cache elements that were created for a given
/// normalization set.
class RooFixedProdPdf : public RooAbsPdf {
public:
   RooFixedProdPdf(RooProdPdf const &prodPdf, RooArgSet const &normSet)
      : RooAbsPdf(prodPdf.GetName(), prodPdf.GetTitle()), _cache(prodPdf.createCacheElem(&normSet, nullptr)),
        _servers("!servers", "List of servers", this), _prodPdf{prodPdf}
   {
      auto &cache = *_cache;

      // The actual servers for a given normalization set depend on whether the
      // cache is rearranged or not. See RooProdPdf::calculateBatch to see
      // which args in the cache are used directly.
      if (cache._isRearranged) {
         _servers.add(*cache._rearrangedNum);
         _servers.add(*cache._rearrangedDen);
      } else {
         for (std::size_t i = 0; i < cache._partList.size(); ++i) {
            _servers.add(cache._partList[i]);
         }
      }
   }
   RooFixedProdPdf(const RooFixedProdPdf &other, const char *name = nullptr)
      : RooAbsPdf(other, name), _servers("!servers", this, other._servers), _prodPdf{other._prodPdf}
   {
   }
   TObject *clone(const char *newname) const override { return new RooFixedProdPdf(*this, newname); }

   bool selfNormalized() const override { return true; }

   inline bool canComputeBatchWithCuda() const override { return true; }

   void computeBatch(cudaStream_t *stream, double *output, size_t nEvents,
                     RooFit::Detail::DataMap const &dataMap) const override
   {
      _prodPdf.calculateBatch(*_cache, stream, output, nEvents, dataMap);
   }

   ExtendMode extendMode() const override { return _prodPdf.extendMode(); }
   double expectedEvents(const RooArgSet *nset) const override { return _prodPdf.expectedEvents(nset); }

   // Analytical Integration handling
   bool forceAnalyticalInt(const RooAbsArg &dep) const override { return _prodPdf.forceAnalyticalInt(dep); }
   Int_t getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &analVars, const RooArgSet *normSet,
                                 const char *rangeName = nullptr) const override
   {
      return _prodPdf.getAnalyticalIntegralWN(allVars, analVars, normSet, rangeName);
   }
   Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &numVars, const char *rangeName = nullptr) const override
   {
      return _prodPdf.getAnalyticalIntegral(allVars, numVars, rangeName);
   }
   double analyticalIntegralWN(Int_t code, const RooArgSet *normSet, const char *rangeName) const override
   {
      return _prodPdf.analyticalIntegralWN(code, normSet, rangeName);
   }
   double analyticalIntegral(Int_t code, const char *rangeName = nullptr) const override
   {
      return _prodPdf.analyticalIntegral(code, rangeName);
   }

private:
   double evaluate() const override { return _prodPdf.calculate(*_cache); }

   std::unique_ptr<RooProdPdf::CacheElem> _cache;
   RooSetProxy _servers;
   RooProdPdf const &_prodPdf;
};

namespace {

using RooFit::Detail::DataKey;
using ServerLists = std::map<DataKey, std::vector<DataKey>>;

class GraphChecker {
public:
   GraphChecker(RooAbsArg const &topNode)
   {
      // Get the list of servers for each node by data key.
      {
         RooArgList nodes;
         topNode.treeNodeServerList(&nodes, nullptr, true, true, false, true);
         RooArgSet nodesSet{nodes};
         for (RooAbsArg *node : nodesSet) {
            _serverLists[node];
            for (RooAbsArg *server : node->servers()) {
               _serverLists[node].push_back(server);
            }
         }
      }
      for (auto &item : _serverLists) {
         auto &l = item.second;
         std::sort(l.begin(), l.end());
         l.erase(std::unique(l.begin(), l.end()), l.end());
      }
   }

   bool dependsOn(DataKey arg, DataKey testArg)
   {

      std::pair<DataKey, DataKey> p{arg, testArg};

      auto found = _results.find(p);
      if (found != _results.end())
         return found->second;

      if (arg == testArg)
         return true;

      auto const &serverList = _serverLists.at(arg);

      // Next test direct dependence
      auto foundServer = std::find(serverList.begin(), serverList.end(), testArg);
      if (foundServer != serverList.end()) {
         _results.emplace(p, true);
         return true;
      }

      // If not, recurse
      for (auto const &server : serverList) {
         bool t = dependsOn(server, testArg);
         _results.emplace(std::pair<DataKey, DataKey>{server, testArg}, t);
         if (t) {
            return true;
         }
      }

      _results.emplace(p, false);
      return false;
   }

private:
   ServerLists _serverLists;
   std::map<std::pair<DataKey, DataKey>, bool> _results;
};

void treeNodeServerListAndNormSets(const RooAbsArg &arg, RooAbsCollection &list, RooArgSet const &normSet,
                                   std::unordered_map<DataKey, RooArgSet *> &normSets, GraphChecker const &checker)
{
   if (normSets.find(&arg) != normSets.end())
      return;

   list.add(arg, true);

   // normalization sets only need to be added for pdfs
   if (dynamic_cast<RooAbsPdf const *>(&arg)) {
      normSets.insert({&arg, new RooArgSet{normSet}});
   }

   // Recurse if current node is derived
   if (arg.isDerived() && !arg.isFundamental()) {
      for (const auto server : arg.servers()) {

         if (!server->isValueServer(arg)) {
            continue;
         }

         auto differentSet = arg.fillNormSetForServer(normSet, *server);
         if (differentSet) {
            differentSet->sort();
         }

         auto &serverNormSet = differentSet ? *differentSet : normSet;

         // Make sure that the server is not already part of the computation
         // graph with a different normalization set.
         auto found = normSets.find(server);
         if (found != normSets.end()) {
            if (found->second->size() != serverNormSet.size() || !serverNormSet.hasSameLayout(*found->second)) {
               std::stringstream ss;
               ss << server->ClassName() << "::" << server->GetName()
                  << " is requested to be evaluated with two different normalization sets in the same model!";
               ss << " This is not supported yet. The conflicting norm sets are:\n    RooArgSet";
               serverNormSet.printValue(ss);
               ss << " requested by " << arg.ClassName() << "::" << arg.GetName() << "\n    RooArgSet";
               found->second->printValue(ss);
               ss << " first requested by other client";
               auto errMsg = ss.str();
               oocoutE(server, Minimization) << errMsg << std::endl;
               throw std::runtime_error(errMsg);
            }
            continue;
         }

         treeNodeServerListAndNormSets(*server, list, serverNormSet, normSets, checker);
      }
   }
}

std::vector<std::unique_ptr<RooAbsArg>> unfoldIntegrals(RooAbsArg const &topNode, RooArgSet const &normSet,
                                                        std::unordered_map<DataKey, RooArgSet *> &normSets,
                                                        RooArgSet &replacedArgs, RooArgSet &newArgs)
{
   std::vector<std::unique_ptr<RooAbsArg>> newNodes;

   // No normalization set: we don't need to create any integrals
   if (normSet.empty())
      return newNodes;

   GraphChecker checker{topNode};

   RooArgSet nodes;
   // The norm sets are sorted to compare them for equality more easliy
   RooArgSet normSetSorted{normSet};
   normSetSorted.sort();
   treeNodeServerListAndNormSets(topNode, nodes, normSetSorted, normSets, checker);

   // Clean normsets of the variables that the arg does not depend on
   for (auto &item : normSets) {
      if (!item.second || item.second->empty())
         continue;
      auto actualNormSet = new RooArgSet{};
      for (auto *narg : *item.second) {
         if (checker.dependsOn(item.first, narg))
            // Add the arg from the actual node list in the computation graph.
            // Like this, we don't accidentally add internal variable clones
            // that the client args returned. Looking this up is fast because
            // of the name pointer hash map optimization.
            actualNormSet->add(*nodes.find(*narg));
      }
      delete item.second;
      item.second = actualNormSet;
   }

   // Function to `oldArg` with `newArg` in the computation graph.
   auto replaceArg = [&](RooAbsArg &newArg, RooAbsArg const &oldArg) {
      const std::string attrib = std::string("ORIGNAME:") + oldArg.GetName();

      newArg.setAttribute(attrib.c_str());

      RooArgList newServerList{newArg};

      RooArgList originalClients;
      for (auto *client : oldArg.clients()) {
         if (nodes.containsInstance(*client)) {
            originalClients.add(*client);
         }
      }
      for (auto *client : originalClients) {
         if (dynamic_cast<RooAbsCachedPdf *>(client))
            continue;
         client->redirectServers(newServerList, false, true);
      }

      replacedArgs.add(oldArg);
      newArgs.add(newArg);

      newArg.setAttribute(attrib.c_str(), false);
   };

   // Replace all pdfs that need to be normalized with a pdf wrapper that
   // applies the right normalization.
   for (RooAbsArg *node : nodes) {
      if (auto pdf = dynamic_cast<RooAbsPdf *>(node)) {
         RooArgSet const &currNormSet = *normSets.at(pdf);

         if (currNormSet.empty())
            continue;

         // The call to getVal() sets up cached states for this normalization
         // set, which is important in case this pdf is also used by clients
         // using the getVal() interface (without this, test 28 in stressRooFit
         // is failing for example).
         pdf->getVal(currNormSet);

         if (pdf->selfNormalized() && !dynamic_cast<RooAbsCachedPdf *>(pdf))
            continue;

         auto normalizedPdf = std::make_unique<RooNormalizedPdf>(*pdf, currNormSet);

         replaceArg(*normalizedPdf, *pdf);

         newNodes.emplace_back(std::move(normalizedPdf));
      }
   }

   for (RooAbsArg *node : nodes) {

      if (auto prodPdf = dynamic_cast<RooProdPdf *>(node)) {
         RooArgSet const &currNormSet = *normSets.at(prodPdf);
         auto normalizedPdf = std::make_unique<RooFixedProdPdf>(*prodPdf, currNormSet);

         replaceArg(*normalizedPdf, *prodPdf);

         newNodes.emplace_back(std::move(normalizedPdf));

         continue;
      }
   }

   return newNodes;
}

std::unique_ptr<RooAbsArg> compile(RooAbsArg &arg, RooArgSet const &normSet, bool skipCompiling = false)
{
   std::unique_ptr<RooAbsArg> newArg;
   auto *pdf = dynamic_cast<RooAbsPdf *>(&arg);
   if (!skipCompiling && pdf && !normSet.empty() && !(pdf->selfNormalized() && !dynamic_cast<RooAbsCachedPdf *>(pdf)) &&
       !arg.getAttribute("_COMPILED")) {
      newArg = std::make_unique<RooNormalizedPdf>(*pdf, normSet);
      arg.setAttribute("_COMPILED"); // to avoid that this arg is compiled again
   } else {
      newArg = std::unique_ptr<RooAbsArg>{static_cast<RooAbsArg *>(arg.Clone())};
      arg.setAttribute("_COMPILED", false); // we are dont with this node and can reset the _COMPILED attribute
   }
   const std::string attrib = std::string("ORIGNAME:") + arg.GetName();
   newArg->setAttribute(attrib.c_str());
   return newArg;
}

void addServerClonesToList(RooAbsArg &arg, RooArgList &clonedArgs,
                           std::unordered_map<TNamed const *, RooAbsArg *> &clonedArgsSet, RooArgSet const &normSet,
                           std::vector<std::unique_ptr<RooArgList>> &working, std::size_t recursionDepth)
{
   const bool skipCompiling = dynamic_cast<RooAbsCachedPdf const *>(&arg);
   if (working.size() < recursionDepth + 1) {
      working.emplace_back(std::make_unique<RooArgList>());
   }
   RooArgList &serverClones = *working[recursionDepth];
   serverClones.clear();
   for (const auto server : arg.servers()) {
      auto existingServerClone = clonedArgsSet.find(server->namePtr());
      if (existingServerClone != clonedArgsSet.end()) {
         serverClones.add(*existingServerClone->second);
      } else if (!server->isFundamental() || normSet.find(*server)) {
         std::unique_ptr<RooAbsArg> serverClone = compile(*server, normSet, skipCompiling);
         clonedArgsSet.emplace(serverClone->namePtr(), serverClone.get());
         serverClones.add(*serverClone);
         addServerClonesToList(*serverClone, clonedArgs, clonedArgsSet, normSet, working, recursionDepth + 1);
         clonedArgs.addOwned(std::move(serverClone));
      }
   }
   arg.redirectServers(serverClones, false, true);
}

} // namespace

namespace RooFit {

namespace Detail {

std::unique_ptr<RooAbsArg> compileForNormSetImpl(RooAbsArg const &arg, RooArgSet const &normSet)
{
   const bool skipCompiling = dynamic_cast<RooAbsCachedPdf const *>(&arg);
   std::unique_ptr<RooAbsArg> head = compile(const_cast<RooAbsArg &>(arg), normSet, skipCompiling);
   RooArgList clonedArgs;
   std::unordered_map<TNamed const *, RooAbsArg *> clonedArgsSet;
   std::vector<std::unique_ptr<RooArgList>> working;
   addServerClonesToList(*head, clonedArgs, clonedArgsSet, normSet, working, 0);

   head->addOwnedComponents(std::move(clonedArgs));

   return head;
}

} // namespace Detail

} // namespace RooFit

/// \class NormalizationIntegralUnfolder
/// \ingroup Roofitcore
///
/// A NormalizationIntegralUnfolder takes the top node of a computation graph
/// and a normalization set for its constructor. The normalization integrals
/// for the PDFs in that graph will be created, and placed into the computation
/// graph itself, rewiring the existing RooAbsArgs.
///
/// Note that for evaluation, the original topNode should not be used anymore,
/// because if it is a pdf there is now a new normalized pdf wrapping it,
/// serving as the new top node. This normalized top node can be retreived by
/// NormalizationIntegralUnfolder::arg().

RooFit::NormalizationIntegralUnfolder::NormalizationIntegralUnfolder(RooAbsArg const &topNode, RooArgSet const &normSet)
   : _topNodeWrapper{std::make_unique<RooAddition>("_dummy", "_dummy", RooArgList{topNode})}, _normSetWasEmpty{
                                                                                                 normSet.empty()}
{
   auto ownedArgs = unfoldIntegrals(*_topNodeWrapper, normSet, _normSets, _replacedArgs, _newArgs);
   for (std::unique_ptr<RooAbsArg> &arg : ownedArgs) {
      _topNodeWrapper->addOwnedComponents(std::move(arg));
   }
   _arg = &static_cast<RooAddition &>(*_topNodeWrapper).list()[0];
}

RooFit::NormalizationIntegralUnfolder::~NormalizationIntegralUnfolder()
{
   for (auto &item : _normSets) {
      delete item.second;
   }
}
