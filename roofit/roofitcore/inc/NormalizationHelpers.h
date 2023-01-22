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
   CompileContext(RooArgSet const &topLevelNormSet);

   ~CompileContext();

   RooAbsArg *compile(RooAbsArg &arg, RooAbsArg &owner, RooArgSet const &normSet);

   void compileServers(RooAbsArg &arg, RooArgSet const &normSet);

private:
   void add(RooAbsArg &arg);
   RooAbsArg *find(RooAbsArg &arg) const;

   RooArgSet const &_topLevelNormSet;
   std::unordered_map<TNamed const *, RooAbsArg *> _clonedArgsSet;
};

template <class T>
std::unique_ptr<T> compileForNormSet(T const &arg, RooArgSet const &normSet)
{
   RooFit::CompileContext ctx{normSet};
   std::unique_ptr<RooAbsArg> head = arg.compileForNormSet(normSet, ctx);
   return std::unique_ptr<T>{static_cast<T *>(head.release())};
}

} // namespace RooFit

#endif
