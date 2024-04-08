#ifndef TMVA_RBDT
#define TMVA_RBDT

#include <Rtypes.h>

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace TMVA {

namespace Experimental {

class RBDT {
public:
   // The floating point number type that will be used to accept features and store cut values
   typedef float FeatureType;
   // Tue floating point number type that the individual trees return their responses in
   typedef float TreeResponseType;
   // The floating point number type that is used to sum the individual tree responses
   typedef float TreeEnsembleResponseType;
   // This integer type stores the indices of the feature employed in each cut.
   typedef unsigned int CutIndexType;

   RBDT() = default;

   /// Construct backends from model in ROOT file
   RBDT(const std::string &key, const std::string &filename);

   inline TreeEnsembleResponseType operator()(const FeatureType *array) const { return evaluateBinary(array); }

   template <int nClasses>
   std::array<TreeEnsembleResponseType, nClasses> softmax(const FeatureType *array) const
   {
      // static softmax interface: no manual memory allocation, but requires to know nClasses at compile time
      static_assert(nClasses >= 3, "nClasses should be >= 3");
      std::array<TreeEnsembleResponseType, nClasses> out{};
      evaluate(array, out.data(), nClasses, baseScore_);
      softmaxTransformInplace(out.data(), nClasses);
      return out;
   }

   // dynamic softmax interface with manually allocated std::vector: simple but inefficient
   std::vector<TreeEnsembleResponseType> softmax(const FeatureType *array) const;

   // softmax interface that is not a pure function, but no manual allocation and no compile-time knowledge needed
   void softmax(const FeatureType *array, TreeEnsembleResponseType *out) const;

   int nClasses() const { return baseResponses_.size() > 2 ? baseResponses_.size() : 2; }

   static RBDT load_txt(std::string const &txtpath, std::vector<std::string> &features, int nClasses = 2);
   static RBDT load_txt(std::istream &is, std::vector<std::string> &features, int nClasses = 2);

   std::vector<int> rootIndices_;
   std::vector<CutIndexType> cutIndices_;
   std::vector<FeatureType> cutValues_;
   std::vector<int> leftIndices_;
   std::vector<int> rightIndices_;
   std::vector<TreeResponseType> responses_;
   std::vector<int> treeNumbers_;
   std::vector<TreeEnsembleResponseType> baseResponses_;
   TreeEnsembleResponseType baseScore_ = 0.0;

private:
   void evaluate(const FeatureType *array, TreeEnsembleResponseType *out, int nOut) const;

   static void softmaxTransformInplace(TreeEnsembleResponseType *out, int nOut);

   TreeEnsembleResponseType evaluateBinary(const FeatureType *array) const;

   ClassDefNV(RBDT, 1);
};

} // namespace Experimental

} // namespace TMVA

#endif // TMVA_RBDT
