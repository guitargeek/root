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
#include <RooMsgService.h>
#include <RooProdPdf.h>

#include "RooFitImplHelpers.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
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
/// The factors of a RooProdPdf that are constraint terms with respect to a set
/// of observables, given by the name pointers of the variables.
///
/// A factor is a constraint term if it depends on none of the observables,
/// can't be extended, and depends on at least one parameter of the factors
/// that are no such candidates. If `parameters` is given, only these
/// parameters count.
std::vector<ConstraintFactor> findConstraintFactors(RooProdPdf const &prod,
                                                    std::unordered_set<TNamed const *> const &observables,
                                                    std::unordered_set<TNamed const *> const *parameters)
{
   std::vector<RooAbsPdf *> candidates;
   std::vector<RooArgSet> candidateVars;
   std::unordered_set<TNamed const *> otherParams;

   for (auto *factor : static_range_cast<RooAbsPdf *>(prod.pdfList())) {
      RooArgSet vars;
      factor->getParameters(nullptr, vars);
      const bool dependsOnObservables =
         std::any_of(vars.begin(), vars.end(), [&](RooAbsArg *v) { return observables.count(v->namePtr()) > 0; });
      if (dependsOnObservables || factor->canBeExtended()) {
         for (RooAbsArg *v : vars) {
            if (!observables.count(v->namePtr()) && (!parameters || parameters->count(v->namePtr()))) {
               otherParams.insert(v->namePtr());
            }
         }
      } else {
         candidates.push_back(factor);
         candidateVars.emplace_back(std::move(vars));
      }
   }

   std::vector<ConstraintFactor> out;
   for (std::size_t i = 0; i < candidates.size(); ++i) {
      RooArgSet shared;
      for (RooAbsArg *v : candidateVars[i]) {
         if (otherParams.count(v->namePtr())) {
            shared.add(*v);
         }
      }
      if (!shared.empty()) {
         out.push_back({candidates[i], std::move(shared)});
      }
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Throw an exception if a RooProdPdf in the model has a factor that looks
/// like a constraint term: one that depends on none of the variables in the
/// dataset, but on floating parameters of the model that the other factors
/// depend on, too.
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

   RooArgSet components;
   pdf.branchNodeServerList(&components);
   for (auto *prod : dynamic_range_cast<RooProdPdf const *>(components)) {
      if (!prod)
         continue;

      for (ConstraintFactor const &factor : findConstraintFactors(*prod, dataVars, &floatingParams)) {
         std::stringstream ss;
         ss << "The factor \"" << factor.pdf->GetName() << "\" of the RooProdPdf \"" << prod->GetName()
            << "\" depends on none of the observables in the dataset, but on the parameters " << factor.parameters
            << " of the other factors. Such a factor cancels in the normalization of the pdf and does NOT constrain "
               "these parameters in the likelihood. Constraint terms are no part of the pdf: pass them with "
               "ModelConfig::SetExternalConstraints() and use ModelConfig::createNLL() or ModelConfig::fitTo(), "
               "or pass them with the ExternalConstraints() command argument. RooFit::splitOffConstraints() "
               "separates the constraint terms of such a pdf from it.";
         oocoutE(&pdf, InputArguments) << ss.str() << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
}

/// \endcond
