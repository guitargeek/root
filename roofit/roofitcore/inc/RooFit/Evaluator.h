/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2021
 *   Emmanouil Michalainas, CERN 2021
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_Evaluator_h
#define RooFit_Evaluator_h

#include <RooAbsReal.h>
#include <RooFit/EvalContext.h>

#include <RConfig.h>

#include <cstdint>
#include <memory>
#include <utility>

class ChangeOperModeRAII;
class RooAbsArg;

namespace RooBatchCompute {
class AbsBufferManager;
}

namespace RooFit {

struct NodeInfo;

class Evaluator {
public:
   Evaluator(const RooAbsReal &absReal, bool useGPU = false);
   ~Evaluator();

   std::span<const double> run();
   void setInput(std::string const &name, std::span<const double> inputArray, bool isOnDevice);
   RooArgSet getParameters() const;
   void print(std::ostream &os);

   void setOffsetMode(RooFit::EvalContext::OffsetMode);

   std::unique_ptr<ChangeOperModeRAII> setOperModes(RooAbsArg::OperMode opMode);

private:
   void processVariable(NodeInfo &nodeInfo);
   void processCategory(NodeInfo &nodeInfo);
   void setClientsDirty(NodeInfo &nodeInfo);
   std::span<const double> getValHeterogeneous();
   void markGPUNodes();
   void assignToGPU(NodeInfo &info);
   void computeCPUNode(const RooAbsArg *node, NodeInfo &info);
   void setOperMode(RooAbsArg *arg, RooAbsArg::OperMode opMode);
   void syncDataTokens();
   void updateOutputSizes();

   std::unique_ptr<RooBatchCompute::AbsBufferManager> _bufferManager;
   RooAbsReal &_topNode;
   const std::uint64_t _dataTokenOwnerId; // for syncDataTokens(), unique per evaluator
   const bool _useGPU = false;
   int _nEvaluations = 0;
   bool _needToUpdateOutputSizes = false;
   RooFit::EvalContext _evalContextCPU;
   RooFit::EvalContext _evalContextCUDA;
   std::vector<NodeInfo> _nodes;                             // the ordered computation graph
   std::unordered_map<TNamed const *, NodeInfo *> _nodesMap; // for quick lookup of nodes
   // Server objects that got de-duplicated by name in the `_nodes` list,
   // together with the node they are an alias of (for syncDataTokens()).
   std::vector<std::pair<RooAbsArg *, NodeInfo const *>> _aliasedServers;
   std::unique_ptr<ChangeOperModeRAII> _operModeChanges;
};

} // end namespace RooFit

#endif
