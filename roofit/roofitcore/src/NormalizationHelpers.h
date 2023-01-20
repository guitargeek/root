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

#include <RooArgSet.h>

#include <memory>

class RooAbsArg;

namespace RooFit {

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
