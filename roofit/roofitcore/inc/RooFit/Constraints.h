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

#ifndef RooFit_Constraints_h
#define RooFit_Constraints_h

#include <RooFit/ModelConfig.h>

#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>

#include <memory>
#include <vector>

class RooWorkspace;

/// \file RooFit/Constraints.h
/// Helpers for the constraint terms of a model.
///
/// Constraint terms are no part of the pdf of a model: they multiply the
/// likelihood once, and not for each event. A RooFit::ModelConfig has them as
/// external constraints, see RooFit::ModelConfig::SetExternalConstraints(),
/// and a likelihood gets them with the RooFit::ExternalConstraints() command
/// argument. The functions here answer the questions that come with this
/// separation: which constraint terms belong to a pdf, which ones a workspace
/// knows, and how to get the values of their global observables.

namespace RooFit {

/// The name of the workspace set with the constraint terms that were
/// registered in the workspace, see addConstraintTerms().
inline constexpr const char *constraintTermsSetName = "constraintTerms";

RooArgSet
constraintsOf(RooAbsArg const &pdf, RooAbsCollection const &constraints, RooArgSet const *observables = nullptr);

std::vector<ModelConfig *> modelConfigs(RooWorkspace const &ws, RooAbsArg const *pdf = nullptr);
RooArgSet constraintTerms(RooWorkspace const &ws);
RooArgSet constraintsOf(RooWorkspace const &ws, RooAbsArg const &pdf);
void addConstraintTerms(RooWorkspace &ws, RooAbsCollection const &constraints);

std::unique_ptr<RooDataSet>
generateGlobalObservables(RooAbsCollection const &constraints, RooArgSet const &globalObservables, int nEvents = 1);

/// The result of splitOffConstraints().
struct SplitPdf {
   /// The pdf without its constraint terms, or a `nullptr` if the pdf has
   /// none. It shares the nodes that didn't change with the original pdf.
   std::unique_ptr<RooAbsPdf> pdf;
   /// The constraint terms, which are objects of the original pdf.
   RooArgSet constraints;
};

SplitPdf splitOffConstraints(RooAbsPdf const &pdf, RooArgSet const &observables);

} // namespace RooFit

#endif
