#ifndef TMVA_RBDT
#define TMVA_RBDT

#include <TFile.h>

//std::unique_ptr<TFile> file{TFile::Open(filename.c_str(),"READ")};
//if (!file || file->IsZombie()) {
   //throw std::runtime_error("Failed to open input file " + filename);
//}
//auto numOutputs = Internal::GetObjectSafe<std::vector<int>>(file.get(), filename, key + "/num_outputs");
//fNumOutputs = numOutputs->at(0);
//delete numOutputs;

//// Get objective and decide whether to normalize output nodes for example in the multiclass case
//auto objective = Internal::GetObjectSafe<std::string>(file.get(), filename, key + "/objective");
//if (objective->compare("softmax") == 0)
   //fNormalizeOutputs = true;
//else
   //fNormalizeOutputs = false;
//delete objective;
//file->Close();

#include <Rtypes.h>

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace fastforest {

// The floating point number type that will be used to accept features and store cut values
typedef float FeatureType;
// Tue floating point number type that the individual trees return their responses in
typedef float TreeResponseType;
// The floating point number type that is used to sum the individual tree responses
typedef float TreeEnsembleResponseType;
// This integer type stores the indices of the feature employed in each cut.
// Set to `unsigned char` for most compact fastforest ofjects if you have less than 256 features.
typedef unsigned int CutIndexType;

// The base response you have to use with older XGBoost versions might be
// zero, so try to explicitely pass zero to the model evaluation if the
// results from this library are incorrect.
const TreeEnsembleResponseType defaultBaseResponse = 0.5;

namespace details {

void softmaxTransformInplace(TreeEnsembleResponseType *out, int nOut);

}

class FastForest {
public:
   inline TreeEnsembleResponseType
   operator()(const FeatureType *array, TreeEnsembleResponseType baseResponse = defaultBaseResponse) const
   {
      return evaluateBinary(array, baseResponse);
   }

   template <int nClasses>
   std::array<TreeEnsembleResponseType, nClasses>
   softmax(const FeatureType *array, TreeEnsembleResponseType baseResponse = defaultBaseResponse) const
   {
      // static softmax interface: no manual memory allocation, but requires to know nClasses at compile time
      static_assert(nClasses >= 3, "nClasses should be >= 3");
      std::array<TreeEnsembleResponseType, nClasses> out{};
      evaluate(array, out.data(), nClasses, baseResponse);
      details::softmaxTransformInplace(out.data(), nClasses);
      return out;
   }

   // dynamic softmax interface with manually allocated std::vector: simple but inefficient
   std::vector<TreeEnsembleResponseType>
   softmax(const FeatureType *array, TreeEnsembleResponseType baseResponse = defaultBaseResponse) const;

   // softmax interface that is not a pure function, but no manual allocation and no compile-time knowledge needed
   void softmax(const FeatureType *array, TreeEnsembleResponseType *out,
                TreeEnsembleResponseType baseResponse = defaultBaseResponse) const;

   int nClasses() const { return baseResponses_.size() > 2 ? baseResponses_.size() : 2; }

   std::vector<int> rootIndices_;
   std::vector<CutIndexType> cutIndices_;
   std::vector<FeatureType> cutValues_;
   std::vector<int> leftIndices_;
   std::vector<int> rightIndices_;
   std::vector<TreeResponseType> responses_;
   std::vector<int> treeNumbers_;
   std::vector<TreeEnsembleResponseType> baseResponses_;

private:
   void evaluate(const FeatureType *array, TreeEnsembleResponseType *out, int nOut,
                 TreeEnsembleResponseType baseResponse) const;

   TreeEnsembleResponseType evaluateBinary(const FeatureType *array, TreeEnsembleResponseType baseResponse) const;

   ClassDefNV(FastForest, 1);
};

FastForest load_txt(std::string const &txtpath, std::vector<std::string> &features, int nClasses = 2);
FastForest load_txt(std::istream &is, std::vector<std::string> &features, int nClasses = 2);

} // namespace fastforest

#endif // TMVA_RBDT
