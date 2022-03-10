/*
 * Project: RooFit
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef roofit_roofit_AsimovDataTools_h
#define roofit_roofit_AsimovDataTools_h

#include <RooDataSet.h>

class RooAbsPdf;
class RooArgSet;

namespace RooFit {

std::unique_ptr<RooDataSet>
generateAsimovDataset(RooAbsPdf const &pdf, RooArgSet const &observables, int printLevel = 1);

}

#endif
