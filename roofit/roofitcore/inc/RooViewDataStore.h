#ifndef RooViewDataStore_h
#define RooViewDataStore_h

#include "RooAbsDataStore.h"
#include "RooVectorDataStore.h"

class RooViewDataStore : public RooAbsDataStore {
public:
   RooViewDataStore() {}
   RooViewDataStore(std::string_view name, std::string_view title, const RooArgSet &vars, int numEntries,
                    std::vector<double *> const &dataReal);
   RooViewDataStore(const RooViewDataStore &other, const char *newname = 0) : RooAbsDataStore(other, newname)
   {
      _isWeighted = other._isWeighted;
      _numEntries = other._numEntries;
      _sumEntries = other._sumEntries;
      _dataReal = other._dataReal;
      _attachedArgs = other._attachedArgs;
      _currentIndex = other._currentIndex;
   }
   RooViewDataStore(const RooViewDataStore &other, const RooArgSet &vars, const char *newname = 0)
      : RooAbsDataStore(other, vars, newname)
   {
      _isWeighted = other._isWeighted;
      _numEntries = other._numEntries;
      _sumEntries = other._sumEntries;
      _dataReal = other._dataReal;
      _attachedArgs = other._attachedArgs;
      _currentIndex = other._currentIndex;
   }

   WRITE_TSTRING_COMPATIBLE_CONSTRUCTOR(RooViewDataStore)

   RooAbsDataStore *clone(const char *newname = 0) const override { return new RooViewDataStore(*this, newname); }
   RooAbsDataStore *clone(const RooArgSet &vars, const char *newname = 0) const override
   {
      return new RooViewDataStore(*this, vars, newname);
   }

   // Retrieve a row
   using RooAbsDataStore::get;
   const RooArgSet *get(Int_t index) const override;
   Double_t weight() const override { return _isWeighted ? _weights[_currentIndex] : 1.0; }

   Double_t weightError(RooAbsData::ErrorType /*etype*/ = RooAbsData::Poisson) const override { return 0.0; }
   void weightError(Double_t & /*lo*/, Double_t & /*hi*/,
                    RooAbsData::ErrorType /*etype*/ = RooAbsData::Poisson) const override
   {
      throw std::runtime_error("not implemented!");
   }

   Bool_t isWeighted() const override { return _isWeighted; }

   RooBatchCompute::RunContext getBatches(std::size_t /*first*/, std::size_t /*len*/) const override;
   RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const override;

   // General & bookkeeping methods
   Int_t numEntries() const override { return _numEntries; }
   Double_t sumEntries() const override { return _sumEntries; }

   // Buffer redirection routines used in inside RooAbsOptTestStatistics
   void attachBuffers(const RooArgSet &extObs) override;
   void resetBuffers() override
   {
      _attachedArgs.clear();
      _attachedArgs.resize(_vars.size());
   }

   void setExternalWeightArray(const Double_t * /*arrayWgt*/, const Double_t * /*arrayWgtErrLo*/,
                               const Double_t * /*arrayWgtErrHi*/, const Double_t * /*arraySumW2*/) override
   {
      throw std::runtime_error("not implemented!");
   }

   /// As a cache, a RooVectorDataStore will be returned.
   RooAbsDataCache *cache() const override
   {
      if (_cache == nullptr) {
         _cache = std::make_unique<RooVectorDataStore>("_mycache", "_mycache", RooArgSet{});
      }
      return _cache.get();
   }

   // all the stuff that is not supported with the RooViewDataStore
   Int_t fill() override
   {
      throw std::runtime_error(
         "Attempting to call RooViewDataStore::fill(), but the RooViewDataStore doesn't support filling.");
      return 0;
   }
   void reset() override { throw std::runtime_error("RooViewDataStore thinks resetting is not useful."); }
   bool changeObservableName(const char * /*from*/, const char * /*to*/) override
   {
      throw std::runtime_error(
         "RooViewDataStore doesn't implement changeObservableName(). Why would you even want to do that?");
      return false;
   }
   RooAbsArg *addColumn(RooAbsArg & /*var*/, Bool_t /*adjustRange*/ = true) override
   {
      throw std::runtime_error("not implemented!");
      return nullptr;
   }
   RooArgSet *addColumns(const RooArgList & /*varList*/) override
   {
      throw std::runtime_error("not implemented!");
      return nullptr;
   }
   RooAbsDataStore *merge(const RooArgSet & /*allvars*/, std::list<RooAbsDataStore *> /*dstoreList*/) override
   {
      throw std::runtime_error("not implemented!");
      return nullptr;
   }
   void append(RooAbsDataStore & /*other*/) override { throw std::runtime_error("not implemented!"); }
   void loadValues(const RooAbsDataStore * /*tds*/, const RooFormulaVar * /*select*/ = 0,
                   const char * /*rangeName*/ = 0, std::size_t /*nStart*/ = 0,
                   std::size_t /*nStop*/ = std::numeric_limits<std::size_t>::max()) override
   {
      throw std::runtime_error("not implemented!");
   }

private:
   bool _isWeighted = false;

   int _numEntries = 0;
   double _sumEntries = 0;

   std::vector<double *> _dataReal;
   std::vector<RooAbsArg *> _attachedArgs;

   double *_weights = nullptr;

   mutable int _currentIndex = 0;
   mutable std::unique_ptr<RooVectorDataStore> _cache = nullptr; //!

   // ClassDefOverride(RooViewDataStore,0)
};

#endif
