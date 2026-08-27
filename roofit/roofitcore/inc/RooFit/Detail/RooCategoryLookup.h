/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2026
 *
 * Copyright (c) 2026, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_Detail_RooCategoryLookup_h
#define RooFit_Detail_RooCategoryLookup_h

#include <RooAbsReal.h>
#include <RooCategoryProxy.h>
#include <RooListProxy.h>

namespace RooFit::Detail {

class RooCategoryLookup : public RooAbsReal {
public:
   RooCategoryLookup() = default;
   RooCategoryLookup(const char *name, const char *title, RooAbsCategory &cat, RooArgList const &values, double scale);
   RooCategoryLookup(const RooCategoryLookup &other, const char *name = nullptr);

   TObject *clone(const char *newname) const override { return new RooCategoryLookup(*this, newname); }

   RooAbsCategory const &cat() const { return *_cat; }
   RooArgList const &values() const { return _values; }
   double scale() const { return _scale; }

protected:
   double evaluate() const override;
   void doEval(RooFit::EvalContext &) const override;

private:
   RooCategoryProxy _cat;
   RooListProxy _values;
   double _scale = 1.0;

   ClassDefOverride(RooCategoryLookup, 1)
};

} // namespace RooFit::Detail

#endif
