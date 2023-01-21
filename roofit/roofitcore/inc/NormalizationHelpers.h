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

#ifndef RooFit_NormalizationHelpers_h
#define RooFit_NormalizationHelpers_h

#include <TNamed.h>

#include <memory>
#include <unordered_map>

class RooAbsArg;
class RooArgSet;

namespace RooFit {

class CompileContext {
public:
   ~CompileContext();

   RooAbsArg *compile(RooAbsArg &arg, RooAbsArg &caller, RooArgSet const &normSet);

   void compileServers(RooAbsArg &arg, RooArgSet const &normSet);

private:
   void add(RooAbsArg &arg);
   RooAbsArg *find(RooAbsArg &arg) const;
   std::unordered_map<TNamed const *, RooAbsArg *> _clonedArgsSet;
};

namespace Detail {

std::unique_ptr<RooAbsArg> compileForNormSetImpl(RooAbsArg const &arg, RooArgSet const &normSet);

}

template <class T>
std::unique_ptr<T> compileForNormSet(T const &arg, RooArgSet const &normSet)
{
   return std::unique_ptr<T>{static_cast<T *>(Detail::compileForNormSetImpl(arg, normSet).release())};
}

} // namespace RooFit

#endif
