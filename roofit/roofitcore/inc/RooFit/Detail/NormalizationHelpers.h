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

#ifndef RooFit_Detail_NormalizationHelpers_h
#define RooFit_Detail_NormalizationHelpers_h

#include <memory>
#include <string>
#include <unordered_map>

class RooAbsArg;
class RooArgSet;

class TNamed;

namespace RooFit {

namespace Detail {

class CompileContext {
public:
   CompileContext(RooArgSet const &topLevelNormSet);

   ~CompileContext();

   template <class T>
   T *compile(T &arg, RooAbsArg &owner, RooArgSet const &normSet)
   {
      return static_cast<T *>(compileImpl(arg, owner, normSet));
   }

   void compileServers(RooAbsArg &arg, RooArgSet const &normSet);
   void compileServer(RooAbsArg &server, RooAbsArg &arg, RooArgSet const &normSet);

   void markAsCompiled(RooAbsArg &arg) const;
   void markSubtreeAsCompiled(RooAbsArg &arg) const;

   // This information is used for the binned likelihood optimization.
   void setLikelihoodMode(bool flag) { _likelihoodMode = flag; }
   bool likelihoodMode() const { return _likelihoodMode; }
   // Whether the likelihood that this computation graph is compiled for is an
   // extended likelihood. Only meaningful if likelihoodMode() is also true.
   void setExtendedMode(bool flag) { _extendedMode = flag; }
   bool extendedMode() const { return _extendedMode; }
   // Whether the likelihood that this computation graph is compiled for uses
   // bin-by-bin offsetting, i.e. the Offset("bin") option of createNLL().
   // Only meaningful if likelihoodMode() is also true.
   void setBinOffsetMode(bool flag) { _binOffsetMode = flag; }
   bool binOffsetMode() const { return _binOffsetMode; }
   // Whether this computation graph is compiled for a chi-squared fit. Like
   // likelihoodMode(), this enables the simultaneous mixture compilation;
   // the two modes are mutually exclusive.
   void setChi2Mode(bool flag) { _chi2Mode = flag; }
   bool chi2Mode() const { return _chi2Mode; }
   void setBinnedLikelihoodMode(bool flag) { _binnedLikelihoodMode = flag; }
   bool binnedLikelihoodMode() const { return _binnedLikelihoodMode; }
   void setBinWidthFuncFlag(bool flag) { _binWidthFuncFlag = flag; }
   bool binWidthFuncFlag() const { return _binWidthFuncFlag; }

private:
   RooAbsArg *compileImpl(RooAbsArg &arg, RooAbsArg &owner, RooArgSet const &normSet);
   void add(RooAbsArg &arg);
   RooAbsArg *find(RooAbsArg &arg) const;
   bool isMarkedAsCompiled(RooAbsArg const &arg) const;

   RooArgSet const &_topLevelNormSet;
   std::unordered_map<TNamed const *, RooAbsArg *> _clonedArgsSet;
   std::unordered_map<RooAbsArg *, RooAbsArg *> _replacements;

   bool _likelihoodMode = false;
   bool _extendedMode = false;
   bool _chi2Mode = false;
   bool _binOffsetMode = false;
   bool _binnedLikelihoodMode = false;
   bool _binWidthFuncFlag = false;
};

template <class T>
std::unique_ptr<T> compileForNormSet(T const &arg, RooArgSet const &normSet)
{
   RooFit::Detail::CompileContext ctx{normSet};
   return std::unique_ptr<T>{static_cast<T *>(arg.compileForNormSet(normSet, ctx).release())};
}

} // namespace Detail

} // namespace RooFit

#endif
