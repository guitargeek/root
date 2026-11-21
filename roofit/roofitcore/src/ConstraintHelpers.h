/// \cond ROOFIT_INTERNAL

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
#include <RooArgSet.h>

#include <unordered_set>
#include <vector>

class RooAbsData;
class RooAbsPdf;
class RooProdPdf;

std::unique_ptr<RooAbsReal> createConstraintTerm(std::string const &name, RooAbsPdf const &pdf, RooAbsData const &data,
                                                 RooArgSet const *constraints, RooArgSet const *globalObservables,
                                                 const char *globalObservablesTag, bool takeGlobalObservablesFromData);

void checkForConstraintsInProducts(RooAbsPdf const &pdf, RooAbsData const &data);

/// A factor of a RooProdPdf that is a constraint term, see findConstraintFactors().
struct ConstraintFactor {
   RooAbsPdf *pdf = nullptr;
   RooArgSet parameters; ///< The parameters that it shares with the factors that are no constraint terms.
};

std::vector<ConstraintFactor> findConstraintFactors(RooProdPdf const &prod,
                                                    std::unordered_set<TNamed const *> const &observables,
                                                    std::unordered_set<TNamed const *> const *parameters = nullptr);

#endif

/// \endcond
