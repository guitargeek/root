// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RACTION
#define ROOT_RACTION

#include "ROOT/RDF/ColumnReaderUtils.hxx"
#include "ROOT/RDF/GraphNode.hxx"
#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/RColumnReaderBase.hxx"
#include "ROOT/RDF/Utils.hxx" // ColumnNames_t, IsInternalColumn
#include "ROOT/RDF/RLoopManager.hxx"
#include "ROOT/RDF/RVariedAction.hxx"

#include <array>
#include <cstddef> // std::size_t
#include <memory>
#include <string>
#include <vector>

namespace ROOT {
namespace Internal {
namespace RDF {

namespace RDFDetail = ROOT::Detail::RDF;
namespace RDFGraphDrawing = ROOT::Internal::RDF::GraphDrawing;

namespace GraphDrawing {
std::shared_ptr<GraphNode> AddDefinesToGraph(std::shared_ptr<GraphNode> node, const RColumnRegister &colRegister,
                                             const std::vector<std::string> &prevNodeDefines,
                                             std::unordered_map<void *, std::shared_ptr<GraphNode>> &visitedMap);
} // namespace GraphDrawing

// clang-format off
/**
 * \class ROOT::Internal::RDF::RAction
 * \ingroup dataframe
 * \brief A RDataFrame node that produces a result
 * \tparam Helper The action helper type, which implements the concrete action logic (e.g. FillHelper, SnapshotHelper)
 * \tparam PrevNode The type of the parent node in the computation graph
 * \tparam ColumnTypes_t A TypeList with the types of the input columns
 *
 */
// clang-format on
template <typename Helper, typename PrevNode, typename ColumnTypes_t = typename Helper::ColumnTypes_t>
class R__CLING_PTRCHECK(off) RAction : public RActionBase {
   using TypeInd_t = std::make_index_sequence<ColumnTypes_t::list_size>;

   Helper fHelper;
   const std::shared_ptr<PrevNode> fPrevNodePtr;
   PrevNode &fPrevNode;
   /// Column readers per slot and per input column
   std::vector<std::array<RColumnReaderBase *, ColumnTypes_t::list_size>> fValues;

   /// The nth flag signals whether the nth input column is a custom column or not.
   std::array<bool, ColumnTypes_t::list_size> fIsDefine;

public:
   RAction(Helper &&h, const ColumnNames_t &columns, std::shared_ptr<PrevNode> pd, const RColumnRegister &colRegister)
      : RActionBase(pd->GetLoopManagerUnchecked(), columns, colRegister, pd->GetVariations()),
        fHelper(std::forward<Helper>(h)), fPrevNodePtr(std::move(pd)), fPrevNode(*fPrevNodePtr), fValues(GetNSlots())
   {
      fLoopManager->Register(this);

      const auto nColumns = columns.size();
      for (auto i = 0u; i < nColumns; ++i)
         fIsDefine[i] = colRegister.IsDefineOrAlias(columns[i]);
   }

   RAction(const RAction &) = delete;
   RAction &operator=(const RAction &) = delete;

   ~RAction() { fLoopManager->Deregister(this); }

   /**
      Retrieve a wrapper to the result of the action that knows how to merge
      with others of the same type.
   */
   std::unique_ptr<RDFDetail::RMergeableValueBase> GetMergeableValue() const final
   {
      return fHelper.GetMergeableValue();
   }

   void Initialize() final { fHelper.Initialize(); }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      RColumnReadersInfo info{RActionBase::GetColumnNames(), RActionBase::GetColRegister(), fIsDefine.data(),
                              *fLoopManager};
      fValues[slot] = GetColumnReaders(slot, r, ColumnTypes_t{}, info);
      fHelper.InitTask(r, slot);
   }

   template <typename ColType>
   auto GetValueChecked(unsigned int slot, std::size_t readerIdx, Long64_t entry) -> ColType &
   {
      if (auto *val = fValues[slot][readerIdx]->template TryGet<ColType>(entry))
         return *val;

      throw std::out_of_range{"RDataFrame: Action (" + fHelper.GetActionName() +
                              ") could not retrieve value for column '" + fColumnNames[readerIdx] + "' for entry " +
                              std::to_string(entry) +
                              ". You can use the DefaultValueFor operation to provide a default value, or "
                              "FilterAvailable/FilterMissing to discard/keep entries with missing values instead."};
   }

   template <typename... ColTypes, std::size_t... S>
   void CallExec(unsigned int slot, Long64_t entry, TypeList<ColTypes...>, std::index_sequence<S...>)
   {
      ROOT::Internal::RDF::CallGuaranteedOrder{[&](auto &&...args) { return fHelper.Exec(slot, args...); },
                                               GetValueChecked<ColTypes>(slot, S, entry)...};
      (void)entry; // avoid unused parameter warning (gcc 12.1)
   }

   void Run(unsigned int slot, Long64_t entry) final
   {
      // check if entry passes all filters
      if (fPrevNode.CheckFilters(slot, entry))
         CallExec(slot, entry, ColumnTypes_t{}, TypeInd_t{});
   }

   void TriggerChildrenCount() final { fPrevNode.IncrChildrenCount(); }

   /// Clean-up operations to be performed at the end of a task.
   void FinalizeSlot(unsigned int slot) final
   {
      fValues[slot].fill(nullptr);
      fHelper.CallFinalizeTask(slot);
   }

   /// Clean-up and finalize the action result (e.g. merging slot-local results).
   /// It invokes the helper's Finalize method.
   void Finalize() final
   {
      fHelper.Finalize();
      SetHasRun();
   }

   std::shared_ptr<RDFGraphDrawing::GraphNode>
   GetGraph(std::unordered_map<void *, std::shared_ptr<RDFGraphDrawing::GraphNode>> &visitedMap) final
   {
      auto prevNode = fPrevNode.GetGraph(visitedMap);
      const auto &prevColumns = prevNode->GetDefinedColumns();

      // Action nodes do not need to go through CreateFilterNode: they are never common nodes between multiple branches
      const auto nodeType = HasRun() ? RDFGraphDrawing::ENodeType::kUsedAction : RDFGraphDrawing::ENodeType::kAction;
      auto thisNode =
         std::make_shared<RDFGraphDrawing::GraphNode>(fHelper.GetActionName(), visitedMap.size(), nodeType);
      visitedMap[(void *)this] = thisNode;

      auto upmostNode = AddDefinesToGraph(thisNode, GetColRegister(), prevColumns, visitedMap);

      thisNode->AddDefinedColumns(GetColRegister().GenerateColumnNames());
      upmostNode->SetPrevNode(prevNode);
      return thisNode;
   }

   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via RResultPtr::RegisterCallback
   void *PartialUpdate(unsigned int slot) final { return fHelper.CallPartialUpdate(slot); }

   std::unique_ptr<RActionBase> MakeVariedAction(std::vector<void *> &&results) final
   {
      auto &&variations = GetVariations();
      auto &&nVariations = variations.size();
      assert(results.size() == nVariations);

      std::vector<Helper> helpers;
      helpers.reserve(nVariations);

      for (decltype(nVariations) i{}; i < nVariations; i++)
         helpers.emplace_back(fHelper.CallMakeNew(results[i], variations[i]));

      return std::unique_ptr<RActionBase>(new RVariedAction<Helper, PrevNode, ColumnTypes_t>{
         std::move(helpers), GetColumnNames(), fPrevNodePtr, GetColRegister()});
   }

   /**
    * \brief Returns a new action with a cloned helper.
    *
    * \param[in] newResult The result to be filled by the new action (needed to clone the helper).
    * \return A unique pointer to the new action.
    */
   std::unique_ptr<RActionBase> CloneAction(void *newResult) final
   {
      return std::make_unique<RAction>(fHelper.CallMakeNew(newResult), GetColumnNames(), fPrevNodePtr,
                                       GetColRegister());
   }

private:
   ROOT::RDF::SampleCallback_t GetSampleCallback() final { return fHelper.GetSampleCallback(); }
};

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif // ROOT_RACTION
