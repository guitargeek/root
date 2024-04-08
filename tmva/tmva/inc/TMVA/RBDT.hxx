/**********************************************************************************
 * Project: ROOT - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 *                                                                                *
 *                                                                                *
 * Description:                                                                   *
 *                                                                                *
 * Authors:                                                                       *
 *      Stefan Wunsch (stefan.wunsch@cern.ch)                                     *
 *      Jonas Rembser (stefan.wunsch@cern.ch)                                     *
 *                                                                                *
 * Copyright (c) 2024:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (see tmva/doc/LICENSE)                                          *
 **********************************************************************************/

#ifndef TMVA_RBDT
#define TMVA_RBDT

#include <Rtypes.h>
#include <TMVA/RTensor.hxx>

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace TMVA {

namespace Experimental {

class RBDT {
public:
   typedef float Value_t;

   RBDT() = default;

   /// Construct backends from model in ROOT file
   RBDT(const std::string &key, const std::string &filename);

   /// Compute model prediction on a single event
   ///
   /// The method is intended to be used with std::vectors-like containers,
   /// for example RVecs.
   template <typename Vector>
   Vector Compute(const Vector &x)
   {
      std::size_t nOut = baseResponses_.size() > 2 ? baseResponses_.size() : 1;
      Vector y(nOut);
      compute(x.data(), y.data());
      return y;
   }

   /// Compute model prediction on a single event
   inline std::vector<Value_t> Compute(std::vector<Value_t> const &x) { return Compute<std::vector<Value_t>>(x); }

   RTensor<Value_t> Compute(RTensor<Value_t> const &x);

   static RBDT load_txt(std::string const &txtpath, std::vector<std::string> &features, int nClasses = 2);
   static RBDT load_txt(std::istream &is, std::vector<std::string> &features, int nClasses = 2);

   std::vector<int> rootIndices_;
   std::vector<unsigned int> cutIndices_;
   std::vector<Value_t> cutValues_;
   std::vector<int> leftIndices_;
   std::vector<int> rightIndices_;
   std::vector<Value_t> responses_;
   std::vector<int> treeNumbers_;
   std::vector<Value_t> baseResponses_;
   Value_t baseScore_ = 0.0;
   bool logistic_ = false;

private:
   void softmax(const Value_t *array, Value_t *out) const;
   void compute(const Value_t *array, Value_t *out) const;
   Value_t evaluateBinary(const Value_t *array) const;

   ClassDefNV(RBDT, 1);
};

} // namespace Experimental

} // namespace TMVA

#endif // TMVA_RBDT
