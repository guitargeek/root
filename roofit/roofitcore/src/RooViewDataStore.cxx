#include "RooViewDataStore.h"

#include "RooAbsReal.h"
#include "RunContext.h"

RooViewDataStore::RooViewDataStore(std::string_view name, std::string_view title, const RooArgSet &vars, int numEntries,
                                   std::vector<double *> const &dataReal)
   : RooAbsDataStore{name, title, vars}, _isWeighted{false}, _numEntries{numEntries},
     _sumEntries{double(dataReal.size())}, _dataReal{dataReal}
{
   // check if all values are in the range supported by the variables
   for (std::size_t iVar = 0; iVar < vars.size(); ++iVar) {
      auto var = static_cast<RooAbsRealLValue *>(vars[iVar]);
      double minAllowedVal = var->getMin();
      double maxAllowedVal = var->getMax();
      double minVal = (maxAllowedVal - minAllowedVal) / 2.;
      double maxVal = (maxAllowedVal - minAllowedVal) / 2.;
      double *data = dataReal[iVar];
      for (std::size_t i = 0; i < numEntries; ++i) {
         minVal = std::min(minVal, data[i]);
         maxVal = std::max(maxVal, data[i]);
      }
      if (minVal < minAllowedVal || maxVal > maxAllowedVal) {
         auto errorMsg =
            std::string("RooViewDataStore(): values out of range of variable \"") + vars[iVar]->GetName() + "\".";
         throw std::runtime_error(errorMsg);
      }
   }
   _attachedArgs.resize(vars.size());
}

const RooArgSet *RooViewDataStore::get(Int_t index) const
{
   _currentIndex = index;

   for (std::size_t iVar = 0; iVar < _vars.size(); ++iVar) {
      if (_attachedArgs[iVar]) {
         static_cast<RooRealVar *>(_attachedArgs[iVar])->setVal(_dataReal[iVar][index]);
      }
      static_cast<RooRealVar *>(_vars[iVar])->setVal(_dataReal[iVar][index]);
   }

   if (_doDirtyProp) {
      // Raise all dirty flags
      for (auto var : _vars) {
         var->setValueDirty(); // This triggers recalculation of all clients
      }
   }

   if (_cache) {
      _cache->get(index);
   }

   return &_vars;
}

/// Retrieve batches for all observables in this data store.
RooBatchCompute::RunContext RooViewDataStore::getBatches(std::size_t first, std::size_t len) const
{
   RooBatchCompute::RunContext evalData;

   for (std::size_t i = 0; i < _vars.size(); ++i) {
      if (!_attachedArgs[i])
         continue;
      RooSpan<const double> span{_dataReal[i] + first, len};
      assert(span.size() == len);
      evalData.spans.emplace(static_cast<RooAbsReal *>(_attachedArgs[i]), std::move(span));
   }

   // we should also return the batches from the cached variables at some point

   return evalData;
}

RooSpan<const double> RooViewDataStore::getWeightBatch(std::size_t first, std::size_t len) const
{
   if (_isWeighted) {
      RooSpan<const double> span{_weights + first, len};
      assert(span.size() == len);
      return span;
   }
   return {};
}

void RooViewDataStore::attachBuffers(const RooArgSet &extObs)
{
   for (std::size_t i = 0; i < _vars.size(); ++i) {
      RooAbsArg *extArg = extObs.find(_vars[i]->GetName());
      if (extArg) {
         _attachedArgs[i] = extArg;
      }
   }
}
