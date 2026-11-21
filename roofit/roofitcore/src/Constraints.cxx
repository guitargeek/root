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

#include <RooFit/Constraints.h>

#include <RooFit/ModelConfig.h>

#include <RooAbsCategoryLValue.h>
#include <RooExtendPdf.h>
#include <RooGlobalFunc.h>
#include <RooProdPdf.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

#include "ConstraintHelpers.h"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace {

RooArgSet const *registeredConstraintTerms(RooWorkspace const &ws)
{
   auto found = ws.sets().find(RooFit::constraintTermsSetName);
   return found != ws.sets().end() ? &found->second : nullptr;
}

// Removes the constraint terms from a pdf, for RooFit::splitOffConstraints().
class ConstraintStripper {
public:
   ConstraintStripper(RooArgSet const &observables, RooArgSet &constraints) : _constraints{constraints}
   {
      for (RooAbsArg const *observable : observables) {
         _observables.insert(observable->namePtr());
      }
   }

   // The pdf without its constraint terms: the pdf itself if it has none, another node of its tree if only that one
   // is left of a product, or a new node, which is kept in newNodes().
   RooAbsPdf *strip(RooAbsPdf &pdf)
   {
      // A pdf can be used several times, e.g. for several channels: it is replaced by the same pdf everywhere.
      auto found = _stripped.find(&pdf);
      if (found != _stripped.end()) {
         return found->second;
      }
      RooAbsPdf *out = &pdf;
      if (auto *sim = dynamic_cast<RooSimultaneous *>(&pdf)) {
         out = strip(*sim);
      } else if (auto *prod = dynamic_cast<RooProdPdf *>(&pdf)) {
         out = strip(*prod);
      } else if (auto *extended = dynamic_cast<RooExtendPdf *>(&pdf)) {
         out = strip(*extended);
      }
      _stripped[&pdf] = out;
      return out;
   }

   std::vector<std::unique_ptr<RooAbsPdf>> &newNodes() { return _newNodes; }

private:
   template <class Pdf_t>
   Pdf_t *keep(std::unique_ptr<Pdf_t> pdf)
   {
      Pdf_t *out = pdf.get();
      _newNodes.push_back(std::move(pdf));
      return out;
   }

   static std::string newName(RooAbsArg const &arg) { return std::string(arg.GetName()) + "_unconstrained"; }

   // The constraint terms themselves, and not the products that group them.
   void addConstraint(RooAbsPdf &pdf)
   {
      if (auto *prod = dynamic_cast<RooProdPdf *>(&pdf)) {
         for (auto *factor : static_range_cast<RooAbsPdf *>(prod->pdfList())) {
            addConstraint(*factor);
         }
      } else {
         _constraints.add(pdf, /*silent=*/true);
      }
   }

   RooAbsPdf *strip(RooSimultaneous &sim)
   {
      std::map<std::string, RooAbsPdf *> channels;
      bool changed = false;
      for (auto const &state : sim.indexCat()) {
         RooAbsPdf *channel = sim.getPdf(state.first);
         if (!channel) {
            continue; // not every state needs to have a pdf
         }
         RooAbsPdf *stripped = strip(*channel);
         changed |= stripped != channel;
         channels[state.first] = stripped;
      }
      if (!changed) {
         return &sim;
      }
      auto &indexCat = const_cast<RooAbsCategoryLValue &>(sim.indexCat());
      return keep(std::make_unique<RooSimultaneous>(newName(sim).c_str(), sim.GetTitle(), channels, indexCat));
   }

   RooAbsPdf *strip(RooProdPdf &prod)
   {
      RooArgSet removed;
      for (ConstraintFactor const &factor : findConstraintFactors(prod, _observables)) {
         addConstraint(*factor.pdf);
         removed.add(*factor.pdf);
      }
      // The other factors can have constraint terms of their own.
      std::unordered_map<RooAbsArg *, RooAbsArg *> replacements;
      std::vector<RooAbsPdf *> remaining;
      for (auto *factor : static_range_cast<RooAbsPdf *>(prod.pdfList())) {
         if (removed.containsInstance(*factor)) {
            continue;
         }
         RooAbsPdf *stripped = strip(*factor);
         if (stripped != factor) {
            replacements[factor] = stripped;
         }
         remaining.push_back(factor);
      }
      if (removed.empty() && replacements.empty()) {
         return &prod;
      }
      // If only one factor is left of the product, that factor is the pdf, unless the product makes it a conditional
      // pdf.
      if (!removed.empty() && remaining.size() == 1) {
         RooArgSet const *nset = prod.findPdfNSet(*remaining[0]);
         if (!nset || nset->empty()) {
            auto found = replacements.find(remaining[0]);
            return found != replacements.end() ? static_cast<RooAbsPdf *>(found->second) : remaining[0];
         }
      }
      // A copy of the product keeps what it says about the factors that stay, e.g. which ones are conditional.
      std::unique_ptr<RooProdPdf> out{static_cast<RooProdPdf *>(prod.clone(newName(prod).c_str()))};
      out->removePdfs(removed);
      if (!replacements.empty()) {
         out->redirectServers(replacements);
      }
      return keep(std::move(out));
   }

   RooAbsPdf *strip(RooExtendPdf &extended)
   {
      auto &pdf = const_cast<RooAbsPdf &>(extended.pdf());
      RooAbsPdf *stripped = strip(pdf);
      if (stripped == &pdf) {
         return &extended;
      }
      auto &n = const_cast<RooAbsReal &>(extended.getN());
      const char *rangeName = extended.getRangeName() ? extended.getRangeName()->GetName() : nullptr;
      return keep(
         std::make_unique<RooExtendPdf>(newName(extended).c_str(), extended.GetTitle(), *stripped, n, rangeName));
   }

   std::unordered_set<TNamed const *> _observables;
   RooArgSet &_constraints;
   std::vector<std::unique_ptr<RooAbsPdf>> _newNodes;
   std::unordered_map<RooAbsPdf *, RooAbsPdf *> _stripped;
};

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// The constraint terms of a pdf: the ones out of `constraints` that depend on
/// at least one parameter of the pdf.
///
/// The parameters are the variables of the pdf that are not in `observables`,
/// see RooAbsArg::getParameters(), whether they are floating or constant: the
/// constraint term of a constant parameter is still part of the model, e.g.
/// for the sampling of the global observables.
///
/// \param[in] pdf The pdf, a component of it like the channel of a
///            RooSimultaneous, or any other object with parameters. If it is
///            one of the constraint terms itself, it is not part of the result.
/// \param[in] constraints The constraint terms to select from, e.g. the
///            external constraints of a ModelConfig.
/// \param[in] observables The observables of the pdf. They can be left out if
///            no constraint term depends on them.
/// \return The selected constraint terms, in the order of `constraints`.

RooArgSet RooFit::constraintsOf(RooAbsArg const &pdf, RooAbsCollection const &constraints, RooArgSet const *observables)
{
   RooArgSet parameters;
   pdf.getParameters(observables, parameters);

   RooArgSet out;
   for (RooAbsArg *constraint : constraints) {
      if (constraint != &pdf && constraint->dependsOn(parameters)) {
         out.add(*constraint, /*silent=*/true);
      }
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// The ModelConfigs of a workspace. If `pdf` is given, only the ones whose pdf
/// is this object or depends on it, e.g. because it is a channel of the
/// RooSimultaneous of the model.

std::vector<RooFit::ModelConfig *> RooFit::modelConfigs(RooWorkspace const &ws, RooAbsArg const *pdf)
{
   std::vector<ModelConfig *> out;
   for (TObject *obj : ws.allGenericObjects()) {
      auto *mc = dynamic_cast<ModelConfig *>(obj);
      if (!mc) {
         continue;
      }
      if (pdf) {
         RooAbsPdf const *model = mc->GetPdf();
         if (!model || (model->namePtr() != pdf->namePtr() && !model->dependsOn(*pdf))) {
            continue;
         }
      }
      out.push_back(mc);
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// All constraint terms that a workspace knows: the external constraints of
/// its ModelConfigs, and the ones registered with addConstraintTerms().

RooArgSet RooFit::constraintTerms(RooWorkspace const &ws)
{
   RooArgSet out;
   for (ModelConfig const *mc : modelConfigs(ws)) {
      if (RooArgSet const *constraints = mc->GetExternalConstraints()) {
         out.add(*constraints, /*silent=*/true);
      }
   }
   if (RooArgSet const *registered = registeredConstraintTerms(ws)) {
      out.add(*registered, /*silent=*/true);
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// The constraint terms of a pdf according to the workspace that it is in.
///
/// These are the external constraints of the ModelConfigs whose pdf is the
/// given one or depends on it, as far as they depend on a parameter of the
/// pdf, see constraintsOf(RooAbsArg const &, RooAbsCollection const &,
/// RooArgSet const *). A ModelConfig is what says which constraint terms its
/// model has, so no other terms are considered for such a pdf.
///
/// For a pdf that is part of no ModelConfig, they are selected in the same way
/// from all constraint terms that the workspace knows, see
/// constraintTerms(): the ones registered with addConstraintTerms(), and the
/// ones of the ModelConfigs of other pdfs with the same parameters.

RooArgSet RooFit::constraintsOf(RooWorkspace const &ws, RooAbsArg const &pdf)
{
   const std::vector<ModelConfig *> models = modelConfigs(ws, &pdf);
   if (models.empty()) {
      return constraintsOf(pdf, constraintTerms(ws));
   }
   RooArgSet out;
   for (ModelConfig const *mc : models) {
      if (RooArgSet const *constraints = mc->GetExternalConstraints()) {
         out.add(constraintsOf(pdf, *constraints, mc->GetObservables()), /*silent=*/true);
      }
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Register constraint terms in a workspace, in its set with the name
/// RooFit::constraintTermsSetName. Terms that are not in the workspace yet are
/// imported.
///
/// A ModelConfig is what says which constraint terms a model has. The
/// registered ones are for pdfs that are part of no ModelConfig, e.g. while a
/// model is being built: constraintsOf(RooWorkspace const &, RooAbsArg const &)
/// finds them, so that they can be passed to a likelihood or to a ModelConfig
/// that is created later. The ModelConfigs of the workspace are not changed.

void RooFit::addConstraintTerms(RooWorkspace &ws, RooAbsCollection const &constraints)
{
   std::string names;
   for (RooAbsArg const *constraint : constraints) {
      if (!ws.arg(constraint->GetName())) {
         ws.import(*constraint, RooFit::RecycleConflictNodes(), RooFit::Silence());
      }
      names += (names.empty() ? "" : ",") + std::string(constraint->GetName());
   }
   if (!names.empty()) {
      ws.extendSet(constraintTermsSetName, names.c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sample values of the global observables from constraint terms, whose
/// product is the pdf of the global observables. The values of the global
/// observables themselves are not changed: to use the sampled values, assign
/// an entry of the returned dataset to the global observables.
///
/// Global observables that no constraint term depends on are not generated.
/// Returns a `nullptr` if there is nothing to generate.

std::unique_ptr<RooDataSet>
RooFit::generateGlobalObservables(RooAbsCollection const &constraints, RooArgSet const &globalObservables, int nEvents)
{
   if (constraints.empty() || globalObservables.empty()) {
      return nullptr;
   }
   RooProdPdf globalObservablesPdf{"globalObservablesPdf", "", constraints};
   std::unique_ptr<RooArgSet> toGenerate{globalObservablesPdf.getObservables(globalObservables)};
   if (toGenerate->empty()) {
      return nullptr;
   }
   return std::unique_ptr<RooDataSet>{globalObservablesPdf.generate(*toGenerate, nEvents)};
}

////////////////////////////////////////////////////////////////////////////////
/// Separate the constraint terms of a pdf from it, for models in which they are
/// factors of the pdf, like the ones that HistFactory built before the
/// constraint terms became external constraints of the ModelConfig.
///
/// Such a pdf can't be used for a likelihood: the factors cancel in the
/// normalization of the pdf, so RooAbsPdf::createNLL() refuses it. This
/// function returns the pdf without these factors and the constraint terms,
/// which can then be passed to ModelConfig::SetExternalConstraints() or to the
/// RooFit::ExternalConstraints() command argument.
///
/// A constraint term is a factor of a RooProdPdf that depends on none of the
/// observables, can't be extended, and depends on at least one parameter of
/// the other factors of the product. They are looked for in the pdf itself,
/// in the channels of a RooSimultaneous, in the factors of a RooProdPdf, and
/// in the pdf of a RooExtendPdf. A RooProdPdf that groups constraint terms is
/// replaced by its factors.
///
/// The original pdf is not changed. The new pdf shares all nodes with it that
/// have no constraint terms, and owns the ones that had to be created, which
/// are named after the original ones with the suffix `_unconstrained`. A
/// RooProdPdf of which only one factor is left is replaced by that factor.
///
/// \param[in] pdf The pdf of the model.
/// \param[in] observables The observables of the model, e.g. the variables of
///            the dataset.
/// \return The pdf without its constraint terms, and the constraint terms. If
///         the pdf has none, the returned pdf is a `nullptr`.

RooFit::SplitPdf RooFit::splitOffConstraints(RooAbsPdf const &pdf, RooArgSet const &observables)
{
   SplitPdf out;
   ConstraintStripper stripper{observables, out.constraints};
   auto &original = const_cast<RooAbsPdf &>(pdf);
   RooAbsPdf *stripped = stripper.strip(original);
   if (stripped == &original) {
      return out;
   }

   auto &newNodes = stripper.newNodes();
   auto found =
      std::find_if(newNodes.begin(), newNodes.end(), [&](auto const &node) { return node.get() == stripped; });
   if (found != newNodes.end()) {
      out.pdf = std::move(*found);
      newNodes.erase(found);
   } else {
      // What is left of the pdf is one of its own nodes.
      out.pdf.reset(static_cast<RooAbsPdf *>(stripped->clone(stripped->GetName())));
   }
   for (auto &node : newNodes) {
      out.pdf->addOwnedComponents(std::move(node));
   }
   return out;
}
