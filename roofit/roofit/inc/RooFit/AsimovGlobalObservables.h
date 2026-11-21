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

#ifndef RooFit_AsimovGlobalObservables_h
#define RooFit_AsimovGlobalObservables_h

class RooAbsCollection;
class RooArgSet;

namespace RooFit {

void setGlobalObservablesToAsimov(RooAbsCollection const &constraints, RooArgSet const &globalObservables);

} // namespace RooFit

#endif
