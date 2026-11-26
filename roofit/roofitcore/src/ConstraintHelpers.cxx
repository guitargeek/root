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

#include "ConstraintHelpers.h"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooConstraintSum.h>
#include <RooFit/Constraints.h>
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooWorkspace.h>

#include "RooFitImplHelpers.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

std::unique_ptr<RooArgSet> getGlobalObservables(RooAbsPdf const &pdf, RooArgSet const &constraints,
                                                RooArgSet const *globalObservables, const char *globalObservablesTag)
{

   if (globalObservables && globalObservablesTag) {
      // error!
      std::string errMsg = "RooAbsPdf::fitTo: GlobalObservables and GlobalObservablesTag options mutually exclusive!";
      oocoutE(&pdf, Minimization) << errMsg << std::endl;
      throw std::invalid_argument(errMsg);
   }
   if (globalObservables) {
      // pass-through of global observables
      return std::make_unique<RooArgSet>(*globalObservables);
   }

   if (globalObservablesTag) {
      oocoutI(&pdf, Minimization) << "User-defined specification of global observables definition with tag named '"
                                  << globalObservablesTag << "'" << std::endl;
   } else {
      // Neither GlobalObservables nor GlobalObservablesTag has been processed -
      // try if a default tag is defined in the head node Check if head not
      // specifies default global observable tag
      if (auto defaultGlobalObservablesTag = pdf.getStringAttribute("DefaultGlobalObservablesTag")) {
         oocoutI(&pdf, Minimization) << "p.d.f. provides built-in specification of global observables definition "
                                     << "with tag named '" << defaultGlobalObservablesTag << "'" << std::endl;
         globalObservablesTag = defaultGlobalObservablesTag;
      }
   }

   if (globalObservablesTag) {
      // The global observables are variables of the constraint terms, which
      // are not part of the pdf.
      RooArgSet allVars;
      for (RooAbsArg *constraint : constraints) {
         RooArgSet constraintVars;
         constraint->getParameters(nullptr, constraintVars);
         allVars.add(constraintVars, /*silent=*/true);
      }
      return std::unique_ptr<RooArgSet>{allVars.selectByAttrib(globalObservablesTag, true)};
   }

   // no global observables specified
   return nullptr;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Create the parameter constraint sum to add to the negative log-likelihood.
/// \return If there are constraints, returns a pointer to the constraint NLL.
///         Returns a `nullptr` if the parameters are unconstrained.
/// \param[in] name Name of the created RooConstraintSum object.
/// \param[in] pdf The PDF model whose parameters should be constrained. The
///            pdf itself never contains constraint terms: they are no part
///            of the probability density of the observables.
/// \param[in] data Dataset used in the fit with the constraint sum. It is
///            used to figure out which are the observables and also to get the
///            global observables definition and values if they are stored in
///            the dataset.
/// \param[in] constraints Set of constraint terms.
/// \param[in] globalObservables The normalization set for the constraint terms.
///            If it is `nullptr`, the set of all parameters of the pdf will
///            be used as the normalization set.
/// \param[in] globalObservablesTag Alternative to define the normalization set
///            for the constraint terms. All variables of the constraint terms
///            that have the attribute with the tag defined by
///            `globalObservablesTag` are used. The `globalObservables` and
///            `globalObservablesTag` parameters are mutually exclusive,
///            meaning at least one of them has to be `nullptr`.
/// \param[in] takeGlobalObservablesFromData If the dataset should be used to automatically
///            define the set of global observables. If this is the case and the
///            set of global observables is still defined manually with the
///            `globalObservables` or `globalObservablesTag` parameters, the
///            values of all global observables that are not stored in the
///            dataset are taken from the model.
std::unique_ptr<RooAbsReal> createConstraintTerm(std::string const &name, RooAbsPdf const &pdf, RooAbsData const &data,
                                                 RooArgSet const *constraints, RooArgSet const *globalObservables,
                                                 const char *globalObservablesTag, bool takeGlobalObservablesFromData)
{
   if (!constraints || constraints->empty()) {
      // no constraints
      return nullptr;
   }

   oocoutI(&pdf, Minimization) << " Including the following constraint terms in minimization: " << *constraints
                               << std::endl;

   // Identify global observables in the model.
   auto glObs = getGlobalObservables(pdf, *constraints, globalObservables, globalObservablesTag);
   RooArgSet pdfParams;
   if (data.getGlobalObservables() && takeGlobalObservablesFromData) {
      if (!glObs) {
         // There were no global observables specified, but there are some in the
         // dataset. We will just take them from the dataset.
         oocoutI(&pdf, Minimization)
            << "The following global observables have been automatically defined according to the dataset "
            << "which also provides their values: " << *data.getGlobalObservables() << std::endl;
         glObs = std::make_unique<RooArgSet>(*data.getGlobalObservables());
      } else {
         // There are global observables specified by the user and also some in
         // the dataset.
         RooArgSet globalsFromDataset;
         data.getGlobalObservables()->selectCommon(*glObs, globalsFromDataset);
         oocoutI(&pdf, Minimization) << "The following global observables have been defined: " << *glObs << ","
                                     << " with the values of " << globalsFromDataset
                                     << " obtained from the dataset and the other values from the model." << std::endl;
      }
   } else if (glObs) {
      oocoutI(&pdf, Minimization)
         << "The following global observables have been defined and their values are taken from the model: " << *glObs
         << std::endl;
      // In this case we don't take global observables from data
      takeGlobalObservablesFromData = false;
   } else {
      pdf.getParameters(data.get(), pdfParams, false);
      oocoutI(&pdf, Minimization)
         << "The global observables are not defined , normalize constraints with respect to the parameters "
         << pdfParams << std::endl;
      takeGlobalObservablesFromData = false;
   }

   return std::make_unique<RooConstraintSum>(name.c_str(), "nllCons", *constraints, glObs ? *glObs : pdfParams,
                                             takeGlobalObservablesFromData);
}

////////////////////////////////////////////////////////////////////////////////
/// The factors of the RooProdPdfs in a model that are constraint terms with
/// respect to a set of observables, given by the name pointers of the
/// variables.
///
/// A factor is a candidate if it depends on none of the observables and can't
/// be extended. A candidate is a constraint term if it depends on a parameter
/// of the rest of the model: of the model without the candidates, or of a
/// candidate that is a constraint term. This is not limited to the other
/// factors of the same product: in a RooSimultaneous, a constraint term can be
/// a factor of a channel that doesn't depend on the parameter, e.g. if the
/// model was built for a statistics tool that collected the constraint terms
/// from all channels. If `parameters` is given, only these parameters count.
///
/// Such factors cancel in the normalized value of the product, so they don't
/// contribute to a likelihood of the observables.
std::vector<ConstraintFactor> findConstraintFactors(RooAbsPdf const &pdf,
                                                    std::unordered_set<TNamed const *> const &observables,
                                                    std::unordered_set<TNamed const *> const *parameters)
{
   auto isParameter = [&](RooAbsArg const &arg) {
      return arg.isFundamental() && !observables.count(arg.namePtr()) &&
             (!parameters || parameters->count(arg.namePtr()));
   };

   // The candidates, by product.
   std::vector<ConstraintFactor> candidates;
   std::unordered_map<RooAbsArg const *, std::unordered_set<RooAbsArg const *>> candidatesOfProduct;
   RooArgSet components;
   pdf.branchNodeServerList(&components);
   for (auto *prod : dynamic_range_cast<RooProdPdf const *>(components)) {
      if (!prod) {
         continue;
      }
      for (auto *factor : static_range_cast<RooAbsPdf *>(prod->pdfList())) {
         RooArgSet vars;
         factor->getParameters(nullptr, vars);
         const bool dependsOnObservables =
            std::any_of(vars.begin(), vars.end(), [&](RooAbsArg *v) { return observables.count(v->namePtr()) > 0; });
         if (!dependsOnObservables && !factor->canBeExtended()) {
            candidatesOfProduct[prod].insert(factor);
            candidates.push_back({prod, factor, std::move(vars)});
         }
      }
   }
   if (candidates.empty()) {
      return {};
   }

   // The parameters of the model without the candidates: the leaves that can be reached through value servers,
   // without going from a product to one of its candidates. A candidate that is reached in another way is part of
   // the model.
   std::unordered_set<TNamed const *> modelParams;
   {
      std::vector<RooAbsArg const *> stack{&pdf};
      std::unordered_set<RooAbsArg const *> visited;
      while (!stack.empty()) {
         RooAbsArg const *node = stack.back();
         stack.pop_back();
         if (!visited.insert(node).second) {
            continue;
         }
         if (node->isFundamental()) {
            if (isParameter(*node)) {
               modelParams.insert(node->namePtr());
            }
            continue;
         }
         auto found = candidatesOfProduct.find(node);
         for (RooAbsArg *server : node->servers()) {
            if (server->isValueServer(*node) && (found == candidatesOfProduct.end() || !found->second.count(server))) {
               stack.push_back(server);
            }
         }
      }
   }

   // A candidate that shares a parameter with the model is a constraint term, and its parameters are parameters of
   // the model then, e.g. for a constraint term of the width of another constraint term.
   std::vector<bool> isConstraint(candidates.size(), false);
   for (bool changed = true; changed;) {
      changed = false;
      for (std::size_t i = 0; i < candidates.size(); ++i) {
         if (isConstraint[i]) {
            continue;
         }
         RooArgSet shared;
         for (RooAbsArg *v : candidates[i].parameters) {
            if (isParameter(*v) && modelParams.count(v->namePtr())) {
               shared.add(*v);
            }
         }
         if (shared.empty()) {
            continue;
         }
         for (RooAbsArg *v : candidates[i].parameters) {
            if (isParameter(*v)) {
               modelParams.insert(v->namePtr());
            }
         }
         candidates[i].parameters.removeAll();
         candidates[i].parameters.add(shared);
         isConstraint[i] = true;
         changed = true;
      }
   }

   std::vector<ConstraintFactor> out;
   for (std::size_t i = 0; i < candidates.size(); ++i) {
      if (isConstraint[i]) {
         out.push_back(std::move(candidates[i]));
      }
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Throw an exception if a RooProdPdf in the model has a factor that looks
/// like a constraint term, see findConstraintFactors(): one that depends on
/// none of the variables in the dataset, but on floating parameters of the
/// rest of the model.
///
/// Such a factor cancels in the normalized value of the product, so it does
/// not contribute to the likelihood. Constraint terms are no part of the pdf:
/// they have to be passed to the likelihood separately. Without this check,
/// the parameters would be silently unconstrained.
void checkForConstraintsInProducts(RooAbsPdf const &pdf, RooAbsData const &data)
{
   std::unordered_set<TNamed const *> dataVars;
   for (RooAbsArg *arg : *data.get()) {
      dataVars.insert(arg->namePtr());
   }

   // Only the floating parameters of the model can be constrained. A variable
   // that the model doesn't depend on as a parameter, e.g. because it is
   // integrated over in a projection of the product, is not constrained.
   std::unordered_set<TNamed const *> floatingParams;
   {
      RooArgSet params;
      pdf.getParameters(data.get(), params);
      for (RooAbsArg *param : params) {
         if (!param->isConstant()) {
            floatingParams.insert(param->namePtr());
         }
      }
   }

   std::vector<ConstraintFactor> factors = findConstraintFactors(pdf, dataVars, &floatingParams);
   if (factors.empty()) {
      return;
   }
   ConstraintFactor const &factor = factors.front();
   std::stringstream ss;
   ss << "The factor \"" << factor.pdf->GetName() << "\" of the RooProdPdf \"" << factor.product->GetName()
      << "\" depends on none of the observables in the dataset, but on the parameters " << factor.parameters
      << " of the model. Such a factor cancels in the normalization of the product and does NOT constrain these "
         "parameters in the likelihood. Constraint terms are no part of the pdf: pass them with "
         "ModelConfig::SetExternalConstraints() and use ModelConfig::createNLL() or ModelConfig::fitTo(), or pass "
         "them with the ExternalConstraints() command argument. RooFit::splitOffConstraints() separates the "
         "constraint terms of such a pdf from it.";
   oocoutE(&pdf, InputArguments) << ss.str() << std::endl;
   throw std::runtime_error(ss.str());
}

////////////////////////////////////////////////////////////////////////////////
/// Throw an exception if a likelihood of the pdf is created without constraint
/// terms although its workspace knows constraint terms for its parameters,
/// see RooFit::constraintsOf(RooWorkspace const &, RooAbsArg const &).
///
/// Constraint terms are no part of the pdf, so such a likelihood would leave
/// the parameters silently unconstrained. This function is called only if no
/// ExternalConstraints() command argument was given: passing one, also an
/// empty one, says which constraint terms the likelihood has.
void checkForForgottenConstraints(RooAbsPdf const &pdf)
{
   RooWorkspace const *ws = pdf.workspace();
   if (!ws) {
      return;
   }
   const RooArgSet constraints = RooFit::constraintsOf(*ws, pdf);
   if (constraints.empty()) {
      return;
   }
   std::stringstream ss;
   ss << "The workspace \"" << ws->GetName() << "\" has constraint terms for parameters of the pdf \"" << pdf.GetName()
      << "\": " << constraints
      << ", but none were passed to the likelihood. Constraint terms are no part of the pdf, so this likelihood "
         "would leave the parameters unconstrained. Use ModelConfig::createNLL() or ModelConfig::fitTo(), or pass the "
         "constraint terms with the ExternalConstraints() command argument, e.g. "
         "ExternalConstraints(RooFit::constraintsOf(ws, pdf)). To create the likelihood without constraint terms on "
         "purpose, pass ExternalConstraints() with an empty set.";
   oocoutE(&pdf, InputArguments) << ss.str() << std::endl;
   throw std::runtime_error(ss.str());
}

/// \endcond
