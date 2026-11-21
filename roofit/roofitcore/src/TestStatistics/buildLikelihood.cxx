/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <RooFit/TestStatistics/buildLikelihood.h>

#include <RooSimultaneous.h>
#include <RooFit/TestStatistics/RooBinnedL.h>
#include <RooFit/TestStatistics/RooUnbinnedL.h>
#include <RooFit/TestStatistics/RooSubsidiaryL.h>
#include <RooFit/TestStatistics/RooSumL.h>
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooRealSumPdf.h>
#include <RooProdPdf.h>
#include <TClass.h>

#include <memory>

namespace RooFit {
/**
 * \brief Namespace for new RooFit test statistic calculation.
 *
 * RooFit::TestStatistics contains a major refactoring of the RooAbsTestStatistic-RooAbsOptTestStatistic-RooNLLVar
 * inheritance tree into:
 *   1. statistics-based classes on the one hand;
 *   2. calculation/evaluation/optimization based classes on the other hand.
 *
 * The likelihood is the central unit on the statistics side. The RooAbsL class is implemented for four kinds of
 * likelihoods: binned, unbinned, "subsidiary" (an optimization for numerical stability that gathers components like
 * global observables) and "sum" (over multiple components of the other types). These classes provide ways to compute
 * their components in parallelizable chunks that can be used by the calculator classes as they see fit.
 *
 * On top of the likelihood classes, we also provide for convenience a likelihood builder `NLLFactory`. This factory
 * analyzes the pdf and automatically constructs the proper likelihood, built up from the available RooAbsL subclasses.
 * Options, like specifying constraint terms or global observables, can be passed using method chaining. The
 * `NLLFactory::build` method finally returns the constructed likelihood as a RooRealL object that can be fit to using
 * RooMinimizer.
 *
 * The calculator "Wrapper" classes are abstract interfaces. These can be implemented for different kinds of algorithms,
 * or with different kinds of optimization "back-ends" in mind. Two fork-based multi-processing implementations based
 * on RooFit::MultiProcess are available, one to calculate the gradient of the likelihood in parallel and one for the
 * likelihood itself. The likelihood can also be calculated serially.
 *
 * The coupling of all these classes to RooMinimizer is made via the MinuitFcnGrad class, which owns the Wrappers that
 * calculate the likelihood components.
 *
 * More extensive documentation is available at
 * https://github.com/root-project/root/blob/master/roofit/doc/developers/test_statistics.md
 */
namespace TestStatistics {

namespace { // private implementation details

/*
 * \brief Build the subsidiary likelihood from the constraint terms
 *
 * \param[in] pdf Raw pointer to the pdf
 * \param[in] data Raw pointer to the dataset
 * \param[in] constraints Set of constraint pdfs. Constraint terms are no
 * part of the pdf, so this is the only way to add them to the likelihood.
 * \param[in] global_observables
 * Observables that have a constant value, independent of the dataset events. They are the normalization set of the
 * constraint terms. If empty, the constraint terms are normalized over the parameters of the pdf.
 * \param[in] global_observables_tag String that can be set as attribute in the variables of the constraint terms
 * to indicate that it is a global observable. Can be used instead of \p global_observables.
 * \return A unique pointer to a RooSubsidiaryL with all constraint terms, or a nullptr if there are none.
 */
std::unique_ptr<RooSubsidiaryL> buildSubsidiaryL(RooAbsPdf *pdf, RooAbsData *data, RooArgSet const &constraints,
                                                 RooArgSet global_observables,
                                                 std::string const &global_observables_tag)
{
   if (constraints.empty()) {
      return nullptr;
   }

   std::string tag = global_observables_tag;
   if (!tag.empty()) {
      oocoutI(nullptr, Minimization) << "User-defined specification of global observables definition with tag named '"
                                     << tag << "'" << std::endl;
   } else if (global_observables.empty()) {
      // neither global_observables nor global_observables_tag was given - try if a default tag is defined in the head
      // node
      if (const char *defGlobObsTag = pdf->getStringAttribute("DefaultGlobalObservablesTag")) {
         oocoutI(nullptr, Minimization)
            << "p.d.f. provides built-in specification of global observables definition with tag named '"
            << defGlobObsTag << "'" << std::endl;
         tag = defGlobObsTag;
      }
   }
   if (!tag.empty()) {
      // The global observables are variables of the constraint terms, which are not part of the pdf.
      global_observables.removeAll();
      RooArgSet allVars;
      for (RooAbsArg *constraint : constraints) {
         RooArgSet constraintVars;
         constraint->getParameters(nullptr, constraintVars);
         allVars.add(constraintVars, /*silent=*/true);
      }
      global_observables.add(*std::unique_ptr<RooArgSet>{allVars.selectByAttrib(tag.c_str(), true)});
   }

   oocoutI(nullptr, Minimization) << " Including the following constraint terms in minimization: " << constraints
                                  << std::endl;
   RooArgSet pdfParams;
   if (!global_observables.empty()) {
      oocoutI(nullptr, Minimization) << "The following global observables have been defined: " << global_observables
                                     << std::endl;
   } else {
      pdf->getParameters(data->get(), pdfParams, false);
   }
   std::string name("likelihood for pdf ");
   name += pdf->GetName();
   return std::make_unique<RooSubsidiaryL>(name, constraints,
                                           (!global_observables.empty()) ? global_observables : pdfParams);
}

/// Get the binned part of a pdf
///
/// \param pdf Raw pointer to the pdf
/// \return A pdf is binned if it has attribute "BinnedLikelihood" and it is a RooRealSumPdf (or derived class).
///         In this case \p pdf itself is returned, otherwise a nullptr.
RooAbsPdf *getBinnedPdf(RooAbsPdf *pdf)
{
   if (pdf->getAttribute("BinnedLikelihood") && pdf->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
      return pdf;
   }
   return nullptr;
}

} // namespace

/*
 * \brief Build a set of likelihood components to build a likelihood from a simultaneous pdf.
 *
 * \return A vector to RooAbsL unique_ptrs that contain all component binned and/or
 * unbinned likelihoods. Note: subsidiary components are not included; use getConstraintsSet and/or
 * buildSubsidiaryLikelihood to add those.
 */
std::vector<std::unique_ptr<RooAbsL>> NLLFactory::getSimultaneousComponents()
{
   auto sim_pdf = dynamic_cast<RooSimultaneous *>(&_pdf);

   // the rest of this function is an adaptation of RooAbsTestStatistic::initSimMode:

   auto &simCat = const_cast<RooAbsCategoryLValue &>(sim_pdf->indexCat());

   // note: this is valid for simultaneous likelihoods, not for other test statistic types (e.g. chi2) for which this
   // should return true.
   bool process_empty_data_sets = RooAbsL::isExtendedHelper(&_pdf, _extended);

   // Note: important not to use cloned dataset here (possible when this code is run in Roo[...]L ctor), use the
   // original one (which is data_ in Roo[...]L ctors, but data here)
   std::vector<std::unique_ptr<RooAbsData>> dsetList{_data.split(*sim_pdf, process_empty_data_sets)};

   // Count number of used states
   std::size_t N_components = 0;

   for (const auto &catState : simCat) {
      // Retrieve the PDF for this simCat state
      RooAbsPdf *component_pdf = sim_pdf->getPdf(catState.first.c_str());
      auto found = std::find_if(dsetList.begin(), dsetList.end(), [&](auto const &item) {
        return catState.first == item->GetName();
      });
      RooAbsData *dset = found != dsetList.end() ? found->get() : nullptr;

      if (component_pdf && dset && (0. != dset->sumEntries() || process_empty_data_sets)) {
         ++N_components;
      }
   }

   // Allocate arrays
   std::vector<std::unique_ptr<RooAbsL>> components;
   components.reserve(N_components);
   //   _gofSplitMode.resize(N_components);  // not used, Hybrid mode only, see below

   // Create array of regular fit contexts, containing subset of data and single fitCat PDF
   std::size_t n = 0;
   for (const auto &catState : simCat) {
      const std::string &catName = catState.first;
      // Retrieve the PDF for this simCat state
      RooAbsPdf *component_pdf = sim_pdf->getPdf(catName.c_str());
      auto found = std::find_if(dsetList.begin(), dsetList.end(), [&](auto const &item) {
        return catName == item->GetName();
      });
      RooAbsData *dset = found != dsetList.end() ? found->get() : nullptr;

      if (component_pdf && dset && (0. != dset->sumEntries() || process_empty_data_sets)) {
         ooccoutI(nullptr, Fitting) << "getSimultaneousComponents: creating slave calculator #" << n << " for state "
                                    << catName << " (" << dset->numEntries() << " dataset entries)" << std::endl;

         if (getBinnedPdf(component_pdf)) {
            components.push_back(std::make_unique<RooBinnedL>(component_pdf, dset));
         } else {
            components.push_back(std::make_unique<RooUnbinnedL>(component_pdf, dset, _extended, _evalBackend));
         }
         components.back()->setSimCount(N_components);

         // Servers may have been redirected between instantiation and (deferred) initialization

         std::unique_ptr<RooArgSet> actualParams{component_pdf->getParameters(dset)};
         RooArgSet params;
         _pdf.getParameters(_data.get(), params);
         RooArgSet selTargetParams;
         params.selectCommon(*actualParams, selTargetParams);

         // Sanity check that the parameters of the component likelihood are among the parameters that the component
         // pdf has in common with the full model. With a non-legacy evaluation backend, the component pdf is
         // additionally "compiled" for the dataset, which can drop factors that cancel in the normalization (see
         // RooProdPdf::compileForNormSet and RooFixedProdPdf). The component likelihood can then legitimately report
         // fewer parameters than the original component pdf, so only containment can be asserted, not set equality.
         assert([&] {
            std::unique_ptr<RooArgSet> componentParams{components.back()->getParameters()};
            for (auto *param : *componentParams) {
               if (!selTargetParams.find(param->GetName())) {
                  return false;
               }
            }
            return true;
         }());

         ++n;
      } else {
         if ((!dset || (0. != dset->sumEntries() && !process_empty_data_sets)) && component_pdf) {
            ooccoutD(nullptr, Fitting) << "getSimultaneousComponents: state " << catName
                                       << " has no data entries, no slave calculator created" << std::endl;
         }
      }
   }
   oocoutI(nullptr, Fitting) << "getSimultaneousComponents: created " << n << " slave calculators." << std::endl;

   return components;
}

/// Create a likelihood builder for a given pdf and dataset.
/// \param[in] pdf Raw pointer to the pdf
/// \param[in] data Raw pointer to the dataset
NLLFactory::NLLFactory(RooAbsPdf &pdf, RooAbsData &data) : _pdf{pdf}, _data{data} {}

/*
 * \brief Build a likelihood from a pdf + dataset, optionally with a subsidiary likelihood component.
 *
 * This function analyzes the pdf and automatically constructs the proper likelihood, built up from the available
 * RooAbsL subclasses. In essence, this can give 8 conceptually different combinations, based on three questions:
 * 1. Is it a simultaneous pdf?
 * 2. Is the pdf binned?
 * 3. Does the pdf have subsidiary terms?
 * If questions 1 and 3 are answered negatively, this function will either return a RooBinnedL or RooUnbinnedL. In all
 * other cases it returns a RooSumL, which will contain RooBinnedL and/or RooUnbinnedL component(s) and possibly a
 * RooSubsidiaryL component with constraint terms.
 *
 * \return A unique pointer to a RooSubsidiaryL that contains all terms in
 * the pdf that can be calculated separately from the other components in the full likelihood.
 */
std::unique_ptr<RooAbsL> NLLFactory::build()
{
   std::unique_ptr<RooAbsL> likelihood;
   std::vector<std::unique_ptr<RooAbsL>> components;

   if (dynamic_cast<RooSimultaneous const *>(&_pdf)) {
      components = getSimultaneousComponents();
   } else if (auto binnedPdf = getBinnedPdf(&_pdf)) {
      likelihood = std::make_unique<RooBinnedL>(binnedPdf, &_data);
   } else { // unbinned
      likelihood = std::make_unique<RooUnbinnedL>(&_pdf, &_data, _extended, _evalBackend);
   }

   auto subsidiary = buildSubsidiaryL(&_pdf, &_data, _externalConstraints, _globalObservables, _globalObservablesTag);
   if (subsidiary) {
      if (likelihood) {
         components.push_back(std::move(likelihood));
      }
      components.push_back(std::move(subsidiary));
   }
   if (!components.empty()) {
      likelihood = std::make_unique<RooSumL>(&_pdf, &_data, std::move(components), _extended);
   }
   return likelihood;
}

/// \param[in] extended Set extended term calculation on, off or use
///            RooAbsL::Extended::Auto to determine automatically based on the
///            pdf whether to activate or not.
NLLFactory &NLLFactory::Extended(RooAbsL::Extended extended)
{
   _extended = extended;
   return *this;
}

/// \param[in] externalConstraints Set of constraint pdfs. Constraint terms
///            are no part of the pdf, so this is the only way to add them to
///            the subsidiary likelihood.
NLLFactory &NLLFactory::ExternalConstraints(const RooArgSet &externalConstraints)
{
   _externalConstraints.add(externalConstraints);
   return *this;
}

/// \param[in] globalObservables Observables that have a constant value,
///            independent of the dataset events. They are the normalization
///            set of the constraint terms.
NLLFactory &NLLFactory::GlobalObservables(const RooArgSet &globalObservables)
{
   _globalObservables.add(globalObservables);
   return *this;
}

/// \param[in] globalObservablesTag String that can be set as attribute in
///            the variables of the constraint terms to indicate that it is a
///            global observable. Can be used instead of
///            NLLFactory::GlobalObservables().
NLLFactory &NLLFactory::GlobalObservablesTag(const char *globalObservablesTag)
{
   _globalObservablesTag = globalObservablesTag;
   return *this;
}

NLLFactory &NLLFactory::EvalBackend(RooFit::EvalBackend evalBackend)
{
   _evalBackend = evalBackend;
   return *this;
}

} // namespace TestStatistics
} // namespace RooFit
