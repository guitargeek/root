/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2021
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_ConstraintHelpers_h
#define RooFit_ConstraintHelpers_h

#include <RooAbsReal.h>

struct ConstraintsInfo {
   RooArgSet constraints;
   RooArgSet normSet;
};

ConstraintsInfo collectConstraints(RooAbsPdf const &pdf, RooAbsData const &data, RooArgSet const *constrainedParameters,
                                   RooArgSet const *externalConstraints, RooArgSet const *globalObservables,
                                   const char *globalObservablesTag, bool takeGlobalObservablesFromData,
                                   bool removeConstraintsFromPdf);

#endif
