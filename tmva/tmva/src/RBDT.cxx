/**********************************************************************************
 * Project: ROOT - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 *                                                                                *
 *                                                                                *
 * Description:                                                                   *
 *                                                                                *
 * Authors:                                                                       *
 *      Jonas Rembser (stefan.wunsch@cern.ch)                                     *
 *                                                                                *
 * Copyright (c) 2024:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (see tmva/doc/LICENSE)                                          *
 **********************************************************************************/

#include <TMVA/RBDT.hxx>

#include <TFile.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

namespace {

namespace detail {

typedef std::map<int, int> IndexMap;

void correctIndices(std::vector<int>::iterator begin, std::vector<int>::iterator end, IndexMap const &nodeIndices,
                    IndexMap const &leafIndices)
{
   for (std::vector<int>::iterator it = begin; it != end; ++it) {
      if (nodeIndices.count(*it)) {
         *it = nodeIndices.at(*it);
      } else if (leafIndices.count(*it)) {
         *it = -leafIndices.at(*it);
      } else {
         throw std::runtime_error("something is wrong in the node structure");
      }
   }
}

template<class Value_t>
void softmaxTransformInplace(Value_t *out, int nOut)
{
   // Do softmax transformation inplace, mimicing exactly the Softmax function
   // in the src/common/math.h source file of xgboost.
   double norm = 0.;
   Value_t wmax = *out;
   int i = 1;
   for (; i < nOut; ++i) {
      wmax = std::max(out[i], wmax);
   }
   i = 0;
   for (; i < nOut; ++i) {
      Value_t &x = out[i];
      x = std::exp(x - wmax);
      norm += x;
   }
   i = 0;
   for (; i < nOut; ++i) {
      out[i] /= static_cast<float>(norm);
   }
}

} // namespace detail

namespace util {

inline bool isInteger(const std::string &s)
{
   if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
      return false;

   char *p;
   strtol(s.c_str(), &p, 10);

   return (*p == 0);
}

template <class NumericType>
struct NumericAfterSubstrOutput {
   explicit NumericAfterSubstrOutput()
   {
      value = 0;
      found = false;
      failed = true;
   }
   NumericType value;
   bool found;
   bool failed;
   std::string rest;
};

template <class NumericType>
inline NumericAfterSubstrOutput<NumericType> numericAfterSubstr(std::string const &str, std::string const &substr)
{
   std::string rest;
   NumericAfterSubstrOutput<NumericType> output;
   output.rest = str;

   std::size_t found = str.find(substr);
   if (found != std::string::npos) {
      output.found = true;
      std::stringstream ss(str.substr(found + substr.size(), str.size() - found + substr.size()));
      ss >> output.value;
      if (!ss.fail()) {
         output.failed = false;
         output.rest = ss.str();
      }
   }
   return output;
}

std::vector<std::string> split(std::string const &strToSplit, char delimeter)
{
   std::stringstream ss(strToSplit);
   std::string item;
   std::vector<std::string> splittedStrings;
   while (std::getline(ss, item, delimeter)) {
      splittedStrings.push_back(item);
   }
   return splittedStrings;
}

bool exists(std::string const &filename)
{
   if (FILE *file = fopen(filename.c_str(), "r")) {
      fclose(file);
      return true;
   } else {
      return false;
   }
}

} // namespace util

} // namespace

using namespace TMVA::Experimental;

/// Compute model prediction on input RTensor
RTensor<TMVA::Experimental::RBDT::Value_t> TMVA::Experimental::RBDT::Compute(RTensor<Value_t> const &x)
{
   std::size_t nOut = baseResponses_.size() > 2 ? baseResponses_.size() : 1;
   const std::size_t rows = x.GetShape()[0];
   const std::size_t cols = x.GetShape()[1];
   RTensor<Value_t> y({rows, nOut}, MemoryLayout::ColumnMajor);
   std::vector<Value_t> xRow(cols);
   std::vector<Value_t> yRow(nOut);
   for (std::size_t iRow = 0; iRow < rows; ++iRow) {
      for (std::size_t iCol = 0; iCol < cols; ++iCol) {
         xRow[iCol] = x({iRow, iCol});
      }
      compute(xRow.data(), yRow.data());
      for (std::size_t iOut = 0; iOut < nOut; ++iOut) {
         y({iRow, iOut}) = yRow[iOut];
      }
   }
   return y;
}

void TMVA::Experimental::RBDT::softmax(const Value_t *array, Value_t *out) const
{
   std::size_t nOut = baseResponses_.size() > 2 ? baseResponses_.size() : 1;
   if (nOut == 1) {
      throw std::runtime_error(
         "Error in RBDT::softmax : binary classification models don't support softmax evaluation. Plase set "
         "the number of classes in the RBDT-creating function if this is a multiclassification model.");
   }

   for (std::size_t i = 0; i < nOut; ++i) {
      out[i] = baseScore_ + baseResponses_[i];
   }

   int iRootIndex = 0;
   for (std::vector<int>::const_iterator indexIter = rootIndices_.begin(); indexIter != rootIndices_.end();
        ++indexIter) {
      int index = *indexIter;
      do {
         int r = rightIndices_[index];
         int l = leftIndices_[index];
         index = array[cutIndices_[index]] < cutValues_[index] ? l : r;
      } while (index > 0);
      out[treeNumbers_[iRootIndex] % nOut] += responses_[-index];
      ++iRootIndex;
   }

   detail::softmaxTransformInplace(out, nOut);
}

void TMVA::Experimental::RBDT::compute(const Value_t *array, Value_t *out) const
{
   std::size_t nOut = baseResponses_.size() > 2 ? baseResponses_.size() : 1;
   if (nOut > 1) {
      softmax(array, out);
   } else {
      out[0] = evaluateBinary(array);
      if (logistic_) {
         out[0] = 1.0 / (1.0 + std::exp(-out[0]));
      }
   }
}

TMVA::Experimental::RBDT::Value_t TMVA::Experimental::RBDT::evaluateBinary(const Value_t *array) const
{
   Value_t out = baseScore_ + baseResponses_[0];

   for (std::vector<int>::const_iterator indexIter = rootIndices_.begin(); indexIter != rootIndices_.end();
        ++indexIter) {
      int index = *indexIter;
      do {
         int r = rightIndices_[index];
         int l = leftIndices_[index];
         index = array[cutIndices_[index]] < cutValues_[index] ? l : r;
      } while (index > 0);
      out += responses_[-index];
   }

   return out;
}

namespace {

void terminateTree(TMVA::Experimental::RBDT &ff, int &nPreviousNodes, int &nPreviousLeaves,
                   detail::IndexMap &nodeIndices, detail::IndexMap &leafIndices, int &treesSkipped)
{
   using detail::correctIndices;
   correctIndices(ff.rightIndices_.begin() + nPreviousNodes, ff.rightIndices_.end(), nodeIndices, leafIndices);
   correctIndices(ff.leftIndices_.begin() + nPreviousNodes, ff.leftIndices_.end(), nodeIndices, leafIndices);

   if (nPreviousNodes != static_cast<int>(ff.cutValues_.size())) {
      ff.treeNumbers_.push_back(ff.rootIndices_.size() + treesSkipped);
      ff.rootIndices_.push_back(nPreviousNodes);
   } else {
      int treeNumbers = ff.rootIndices_.size() + treesSkipped;
      ++treesSkipped;
      ff.baseResponses_[treeNumbers % ff.baseResponses_.size()] += ff.responses_.back();
      ff.responses_.pop_back();
   }

   nodeIndices.clear();
   leafIndices.clear();
   nPreviousNodes = ff.cutValues_.size();
   nPreviousLeaves = ff.responses_.size();
}

} // namespace

RBDT TMVA::Experimental::RBDT::load_txt(std::string const &txtpath, std::vector<std::string> &features, int nClasses)
{
   const std::string info = "constructing RBDT from " + txtpath + ": ";

   if (!util::exists(txtpath)) {
      throw std::runtime_error(info + "file does not exists");
   }

   std::ifstream file(txtpath.c_str());
   return load_txt(file, features, nClasses);
}

RBDT TMVA::Experimental::RBDT::load_txt(std::istream &file, std::vector<std::string> &features, int nClasses)
{
   const std::string info = "constructing RBDT from istream: ";

   RBDT ff;
   ff.baseResponses_.resize(nClasses <= 2 ? 1 : nClasses);

   int treesSkipped = 0;

   int nVariables = 0;
   std::map<std::string, int> varIndices;
   bool fixFeatures = false;

   if (!features.empty()) {
      fixFeatures = true;
      nVariables = features.size();
      for (int i = 0; i < nVariables; ++i) {
         varIndices[features[i]] = i;
      }
   }

   std::string line;

   detail::IndexMap nodeIndices;
   detail::IndexMap leafIndices;

   int nPreviousNodes = 0;
   int nPreviousLeaves = 0;

   while (std::getline(file, line)) {
      std::size_t foundBegin = line.find("[");
      std::size_t foundEnd = line.find("]");
      if (foundBegin != std::string::npos) {
         std::string subline = line.substr(foundBegin + 1, foundEnd - foundBegin - 1);
         if (util::isInteger(subline) && !ff.responses_.empty()) {
            terminateTree(ff, nPreviousNodes, nPreviousLeaves, nodeIndices, leafIndices, treesSkipped);
         } else if (!util::isInteger(subline)) {
            std::stringstream ss(line);
            int index;
            ss >> index;
            line = ss.str();

            std::vector<std::string> splitstring = util::split(subline, '<');
            std::string const &varName = splitstring[0];
            Value_t cutValue;
            {
               std::stringstream ss1(splitstring[1]);
               ss1 >> cutValue;
            }
            if (!varIndices.count(varName)) {
               if (fixFeatures) {
                  throw std::runtime_error(info + "feature " + varName + " not in list of features");
               }
               varIndices[varName] = nVariables;
               features.push_back(varName);
               ++nVariables;
            }
            int yes;
            int no;
            util::NumericAfterSubstrOutput<int> output = util::numericAfterSubstr<int>(line, "yes=");
            if (!output.failed) {
               yes = output.value;
            } else {
               throw std::runtime_error(info + "problem while parsing the text dump");
            }
            output = util::numericAfterSubstr<int>(output.rest, "no=");
            if (!output.failed) {
               no = output.value;
            } else {
               throw std::runtime_error(info + "problem while parsing the text dump");
            }

            ff.cutValues_.push_back(cutValue);
            ff.cutIndices_.push_back(varIndices[varName]);
            ff.leftIndices_.push_back(yes);
            ff.rightIndices_.push_back(no);
            std::size_t nNodeIndices = nodeIndices.size();
            nodeIndices[index] = nNodeIndices + nPreviousNodes;
         }

      } else {
         util::NumericAfterSubstrOutput<Value_t> output = util::numericAfterSubstr<Value_t>(line, "leaf=");
         if (output.found) {
            std::stringstream ss(line);
            int index;
            ss >> index;
            line = ss.str();

            ff.responses_.push_back(output.value);
            std::size_t nLeafIndices = leafIndices.size();
            leafIndices[index] = nLeafIndices + nPreviousLeaves;
         }
      }
   }
   terminateTree(ff, nPreviousNodes, nPreviousLeaves, nodeIndices, leafIndices, treesSkipped);

   if (nClasses > 2 && (ff.rootIndices_.size() + treesSkipped) % nClasses != 0) {
      std::stringstream ss;
      ss << "Error in RBDT construction : Forest has " << ff.rootIndices_.size()
         << " trees, which is not compatible with " << nClasses << "classes!";
      throw std::runtime_error(ss.str());
   }

   return ff;
}

TMVA::Experimental::RBDT::RBDT(const std::string &key, const std::string &filename)
{
   std::unique_ptr<TFile> file{TFile::Open(filename.c_str(), "READ")};
   if (!file || file->IsZombie()) {
      throw std::runtime_error("Failed to open input file " + filename);
   }
   auto *fromFile = file->Get<TMVA::Experimental::RBDT>(key.c_str());
   if (!fromFile) {
      throw std::runtime_error("No RBDT with name " + key);
   }
   *this = *fromFile;
}
