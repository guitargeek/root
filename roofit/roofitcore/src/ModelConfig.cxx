/*
 * Project: RooFit
 * Authors:
 *   Kyle Cranmer,
 *   Lorenzo Moneta,
 *   Gregory Schott,
 *   Wouter Verkerke,
 *   Sven Kreiss
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/** \class RooStats::ModelConfig
    \ingroup Roostats

ModelConfig is a simple class that holds configuration information specifying how a model
should be used in the context of various RooStats tools.  A single model can be used
in different ways, and this class should carry all that is needed to specify how it should be used.
ModelConfig requires a workspace to be set.

A ModelConfig holds sets of parameters of the likelihood function that have different interpretations:
- **Parameter of interest** Parameters that are measured (*i.e.* fitted).
- **Nuisance parameters** Parameters that are fitted, but their post-fit value is not interesting. Often,
they might be constrained because external knowledge about them exists, *e.g.* from external measurements.
- **Constraint parameters** No direct use in RooFit/RooStats. Can be used by the user for bookkeeping.
- **Observables** Parameters that have been measured externally, *i.e.* they exist in a dataset. These are not fitted,
but read during fitting from the entries of a dataset.
- **Conditional observables** Observables that are not integrated when the normalisation of the PDF is calculated.
See *e.g.* `rf306_condpereventerrors` in the RooFit tutorials.
- **Global observables** Observables that to the fit look like "constant" values, *i.e.* they are not being
fitted and they are not loaded from a dataset, but some knowledge exists that allows to set them to a
specific value. Examples:
-- A signal efficiency measured in a Monte Carlo study.
-- When constraining a parameter \f$ b \f$, the target value (\f$ b_0 \f$) that this parameter is constrained to:
\f[
  \mathrm{Constraint}_b = \mathrm{Gauss}(b_0 \, | \, b, 0.2)
\f]
- **External constraints** The constraint terms of the model, which are multiplied with the likelihood of the
observables. The ModelConfig is the only place that knows about them: the pdf is the probability density of the
observables and contains no constraint terms. The product of the constraint terms is the probability density of the
global observables, which the RooStats tools use to sample them for toy experiments.

To make sure the constraint terms are considered, the likelihood has to be created with ModelConfig::createNLL() or
ModelConfig::fitTo(), and not with the functions of the pdf.

The functions in RooFit/Constraints.h help with the constraint terms: RooFit::constraintsOf() selects the ones of a
component of the model, like a channel, RooFit::generateGlobalObservables() samples the global observables, and
RooFit::splitOffConstraints() converts a pdf that has its constraint terms as factors.
RooFit::setGlobalObservablesToAsimov() in RooFit/AsimovGlobalObservables.h sets the global observables to their
expected values.
*/

#include <RooFit/ModelConfig.h>

#include <RooFit/Constraints.h>

#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace {

void removeConstantParameters(RooAbsCollection &coll)
{
   RooArgSet constSet;
   for (auto const *myarg : static_range_cast<RooRealVar *>(coll)) {
      if (myarg->isConstant())
         constSet.add(*myarg);
   }
   coll.remove(constSet);
}

} // namespace

using std::ostream;

namespace RooStats {

////////////////////////////////////////////////////////////////////////////////
/// Makes sensible guesses of observables, parameters of interest
/// and nuisance parameters if one or multiple have been set by the creator of this ModelConfig.
///
/// Defaults:
/// - Observables: determined from data,
/// - Global observables: explicit obs  -  obs from data  -  constant observables
/// - Parameters of interest: empty,
/// - Nuisance parameters: all parameters except parameters of interest
///
/// We use nullptr to mean not set, so we don't want to fill
/// with empty RooArgSets.

void ModelConfig::GuessObsAndNuisance(const RooArgSet &obsSet, bool printModelConfig)
{

   // observables
   if (!GetObservables()) {
      SetObservables(*std::unique_ptr<RooArgSet>{GetPdf()->getObservables(obsSet)});
   }
   // global observables
   if (!GetGlobalObservables()) {
      RooArgSet co(*GetObservables());
      co.remove(*std::unique_ptr<RooArgSet>{GetPdf()->getObservables(obsSet)});
      removeConstantParameters(co);
      if (!co.empty())
         SetGlobalObservables(co);

      // TODO BUG This does not work as observables with the same name are already in the workspace.
      /*
      RooArgSet o(*GetObservables());
      o.remove(co);
      SetObservables(o);
      */
   }

   // parameters
   //   if (!GetParametersOfInterest()) {
   //      SetParametersOfInterest(RooArgSet());
   //   }
   if (!GetNuisanceParameters()) {
      RooArgSet params;
      GetPdf()->getParameters(&obsSet, params);
      RooArgSet p(params);
      p.remove(*GetParametersOfInterest());
      removeConstantParameters(p);
      if (!p.empty())
         SetNuisanceParameters(p);
   }

   // print Modelconfig as an info message

   if (printModelConfig) {
      std::ostream &oldstream = RooPrintable::defaultPrintStream(&ccoutI(InputArguments));
      Print();
      RooPrintable::defaultPrintStream(&oldstream);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print contents of Model on the default print stream
/// It can be changed using RooPrintable

void ModelConfig::Print(Option_t *) const
{
   ostream &os = RooPrintable::defaultPrintStream();

   os << std::endl << "=== Using the following for " << GetName() << " ===" << std::endl;

   // args
   if (GetObservables()) {
      os << "Observables:             ";
      GetObservables()->Print("");
   }
   if (GetParametersOfInterest()) {
      os << "Parameters of Interest:  ";
      GetParametersOfInterest()->Print("");
   }
   if (GetNuisanceParameters()) {
      os << "Nuisance Parameters:     ";
      GetNuisanceParameters()->Print("");
   }
   if (GetGlobalObservables()) {
      os << "Global Observables:      ";
      GetGlobalObservables()->Print("");
   }
   if (GetConstraintParameters()) {
      os << "Constraint Parameters:   ";
      GetConstraintParameters()->Print("");
   }
   if (GetExternalConstraints()) {
      os << "External Constraints:    ";
      GetExternalConstraints()->Print("");
   }
   if (GetConditionalObservables()) {
      os << "Conditional Observables: ";
      GetConditionalObservables()->Print("");
   }
   if (GetProtoData()) {
      os << "Proto Data:              ";
      GetProtoData()->Print("");
   }

   // pdfs
   if (GetPdf()) {
      os << "PDF:                     ";
      GetPdf()->Print("");
   }
   if (GetPriorPdf()) {
      os << "Prior PDF:               ";
      GetPriorPdf()->Print("");
   }

   // snapshot
   const RooArgSet *snapshot = GetSnapshot();
   if (snapshot) {
      os << "Snapshot:                " << std::endl;
      snapshot->Print("v");
      delete snapshot;
   }

   os << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Sample values of the global observables from the constraint terms of the
/// model, whose product is the pdf of the global observables. The values of
/// the global observables in the model are not changed: to use the sampled
/// values, assign an entry of the returned dataset to the global observables.
///
/// Global observables that no constraint term depends on are not generated.
/// Returns a `nullptr` if the model has no global observables or no constraint
/// terms.

std::unique_ptr<RooDataSet> ModelConfig::generateGlobalObservables(int nEvents) const
{
   const RooArgSet *globalObservables = GetGlobalObservables();
   const RooArgSet *constraints = GetExternalConstraints();
   if (!globalObservables || !constraints) {
      return nullptr;
   }
   return RooFit::generateGlobalObservables(*constraints, *globalObservables, nEvents);
}

namespace {

// The problems of the constraint terms of a ModelConfig, see ModelConfig::validate().
void collectConstraintProblems(RooStats::ModelConfig const &mc, std::vector<std::string> &errors,
                               std::vector<std::string> &warnings)
{
   RooArgSet const *constraints = mc.GetExternalConstraints();
   RooAbsPdf const *pdf = mc.GetPdf();
   if (!constraints || !pdf) {
      return;
   }
   RooArgSet const *observables = mc.GetObservables();
   RooArgSet parameters;
   pdf->getParameters(observables, parameters);
   const std::string prefix = "ModelConfig \"" + std::string(mc.GetName()) + "\": the external constraint \"";
   for (RooAbsArg const *term : *constraints) {
      if (!dynamic_cast<RooAbsPdf const *>(term)) {
         errors.push_back(prefix + term->GetName() + "\" is no pdf.");
      } else if (observables && term->dependsOn(*observables)) {
         RooArgSet dependents;
         term->getObservables(observables, dependents);
         std::stringstream ss;
         ss << prefix << term->GetName() << "\" depends on the observables " << dependents
            << ". A constraint term is a probability density of global observables, which have one value in the "
               "whole likelihood, so it multiplies the likelihood only once. A probability density of the observables "
               "is part of the pdf.";
         errors.push_back(ss.str());
      } else if (!term->dependsOn(parameters)) {
         warnings.push_back(prefix + term->GetName() + "\" depends on none of the parameters of the pdf \"" +
                            pdf->GetName() +
                            "\", so it constrains nothing in this model. Was the model changed after the constraint "
                            "terms were set? Its parameters are fitted only to the constraint term.");
      }
   }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Check that the constraint terms fit the model. The ModelConfig is the only
/// link between the pdf and its constraint terms, so they can get out of sync,
/// e.g. if the pdf is changed in the workspace after the constraint terms were
/// set. These problems are reported:
///
/// - **Errors**: an external constraint that is no pdf, or that depends on an
///   observable. A constraint term is a probability density of global
///   observables, which multiplies the likelihood only once.
/// - **Warnings**: an external constraint that depends on none of the
///   parameters of the pdf. It constrains nothing in this model, and its
///   parameters are fitted only to the constraint term.
///
/// ModelConfig::createNLL() and ModelConfig::fitTo() do the same check: they
/// throw an exception for errors, and report warnings once.
///
/// \return `true` if no problems were found.

bool ModelConfig::validate() const
{
   std::vector<std::string> errors;
   std::vector<std::string> warnings;
   collectConstraintProblems(*this, errors, warnings);
   for (std::string const &error : errors) {
      oocoutE(this, InputArguments) << error << std::endl;
   }
   for (std::string const &warning : warnings) {
      oocoutW(this, InputArguments) << warning << std::endl;
   }
   return errors.empty() && warnings.empty();
}

////////////////////////////////////////////////////////////////////////////////
/// The check of validate() for createNLL() and fitTo(): throws an exception
/// for errors, and reports warnings only the first time.

void ModelConfig::checkConstraintTerms() const
{
   std::vector<std::string> errors;
   std::vector<std::string> warnings;
   collectConstraintProblems(*this, errors, warnings);
   if (!errors.empty()) {
      for (std::string const &error : errors) {
         oocoutE(this, InputArguments) << error << std::endl;
      }
      throw std::runtime_error(errors.front());
   }
   if (!fConstraintWarningsShown) {
      for (std::string const &warning : warnings) {
         oocoutW(this, InputArguments) << warning << std::endl;
      }
      fConstraintWarningsShown = true;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get all variables of the model: the ones of the pdf and the ones of the
/// constraint terms, which are no part of the pdf. This includes the
/// observables and the global observables.

std::unique_ptr<RooArgSet> ModelConfig::getVariables() const
{
   auto vars = std::make_unique<RooArgSet>();
   if (RooAbsPdf *pdf = GetPdf()) {
      pdf->getParameters(nullptr, *vars);
   }
   if (const RooArgSet *constraints = GetExternalConstraints()) {
      for (RooAbsArg *constraint : *constraints) {
         RooArgSet constraintVars;
         constraint->getParameters(nullptr, constraintVars);
         vars->add(constraintVars, /*silent=*/true);
      }
   }
   return vars;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all parameters of the model: the variables of the pdf and of the
/// constraint terms (see getVariables()) that are not in the set of
/// observables, where the variables are matched by name. Like for
/// RooAbsArg::getParameters(), the global observables are among the parameters
/// if they are no observables in the set.

std::unique_ptr<RooArgSet> ModelConfig::getParameters(const RooArgSet &observables) const
{
   std::unique_ptr<RooArgSet> params = getVariables();
   params->remove(observables, /*silent=*/true, /*matchByNameOnly=*/true);
   return params;
}

////////////////////////////////////////////////////////////////////////////////
/// If a workspace already exists in this ModelConfig, this function will do nothing.

void ModelConfig::SetWS(RooWorkspace &ws)
{
   if (!fRefWS) {
      fRefWS = &ws;
   } else {
      coutE(ObjectHandling) << "ModelConfig::SetWS(): workspace already set, not doing anything" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get workspace.

RooWorkspace *ModelConfig::GetWS() const
{
   if (!fRefWS) {
      coutE(ObjectHandling) << "workspace not set" << std::endl;
      return nullptr;
   }
   return fRefWS;
}

////////////////////////////////////////////////////////////////////////////////
/// save snapshot in the workspace
/// and use values passed with the set

void ModelConfig::SetSnapshot(const RooArgSet &set)
{
   if (!GetWS())
      return;

   fSnapshotName = GetName();
   if (!fSnapshotName.empty())
      fSnapshotName += "_";
   fSnapshotName += set.GetName();
   if (!fSnapshotName.empty())
      fSnapshotName += "_";
   fSnapshotName += "snapshot";
   GetWS()->saveSnapshot(fSnapshotName, set, true); // import also the given parameter values
   DefineSetInWS(fSnapshotName.c_str(), set);
}

////////////////////////////////////////////////////////////////////////////////
/// Load the snapshot from ws and return the corresponding set with the snapshot values.
/// User must delete returned RooArgSet.

const RooArgSet *ModelConfig::GetSnapshot() const
{
   if (!GetWS())
      return nullptr;
   if (fSnapshotName.empty())
      return nullptr;
   // calling loadSnapshot will also copy the current parameter values in the workspaces
   // since we do not want to change the model parameters - we restore the previous ones
   if (!GetWS()->set(fSnapshotName))
      return nullptr;
   RooArgSet snapshotVars(*GetWS()->set(fSnapshotName));
   if (snapshotVars.empty())
      return nullptr;
   // make my snapshot which will contain a copy of the snapshot variables
   RooArgSet tempSnapshot;
   snapshotVars.snapshot(tempSnapshot);
   // load snapshot value from the workspace
   if (!(GetWS()->loadSnapshot(fSnapshotName.c_str())))
      return nullptr;
   // by doing this snapshotVars will have the snapshot values - make the snapshot to return
   const RooArgSet *modelSnapshot = dynamic_cast<const RooArgSet *>(snapshotVars.snapshot());
   // restore now the variables of snapshot in ws to their original values
   // need to const cast since assign is not const (but in reality in just assign values and does not change the set)
   // and anyway the set is const
   snapshotVars.assignFast(tempSnapshot);
   return modelSnapshot;
}

////////////////////////////////////////////////////////////////////////////////
/// load the snapshot from ws if it exists

void ModelConfig::LoadSnapshot() const
{
   if (!GetWS())
      return;
   GetWS()->loadSnapshot(fSnapshotName.c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// helper functions to avoid code duplication

void ModelConfig::DefineSetInWS(const char *name, const RooArgSet &set)
{
   if (!GetWS())
      return;

   const RooArgSet *prevSet = GetWS()->set(name);
   if (prevSet) {
      // be careful not to remove passed set in case it is the same updated
      if (prevSet != &set)
         GetWS()->removeSet(name);
   }

   // suppress warning when we re-define a previously defined set (when set == prevSet )
   // and set is not removed in that case
   RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

   GetWS()->defineSet(name, set, true);

   RooMsgService::instance().setGlobalKillBelow(level);
}

////////////////////////////////////////////////////////////////////////////////
/// internal function to import Pdf in WS

void ModelConfig::ImportPdfInWS(const RooAbsPdf &pdf)
{
   if (!GetWS())
      return;

   if (!GetWS()->pdf(pdf.GetName())) {
      RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      GetWS()->import(pdf, RooFit::RecycleConflictNodes());
      RooMsgService::instance().setGlobalKillBelow(level);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// internal function to import data in WS

void ModelConfig::ImportDataInWS(RooAbsData &data)
{
   if (!GetWS())
      return;

   if (!GetWS()->data(data.GetName())) {
      RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      GetWS()->import(data);
      RooMsgService::instance().setGlobalKillBelow(level);
   }
}

////////////////////////////////////////////////////////////////////////////////

bool ModelConfig::SetHasOnlyParameters(const RooArgSet &set, const char *errorMsgPrefix)
{

   RooArgSet nonparams;
   for (auto const *arg : set) {
      if (!arg->isFundamental()) {
         nonparams.add(*arg);
      }
   }

   if (errorMsgPrefix && !nonparams.empty()) {
      std::cout << errorMsgPrefix << " ERROR: specified set contains non-parameters: " << nonparams << std::endl;
   }
   return (nonparams.empty());
}

/// Specify the external constraints.
void ModelConfig::SetExternalConstraints(const RooArgSet &set)
{
   fExtConstraintsName = std::string(GetName()) + "_ExternalConstraints";
   DefineSetInWS(fExtConstraintsName.c_str(), set);
   fConstraintWarningsShown = false;
}

/// Specify the conditional observables.
void ModelConfig::SetConditionalObservables(const RooArgSet &set)
{
   if (!SetHasOnlyParameters(set, "ModelConfig::SetConditionalObservables"))
      return;
   fConditionalObsName = std::string(GetName()) + "_ConditionalObservables";
   DefineSetInWS(fConditionalObsName.c_str(), set);
}

/// Specify the global observables.
void ModelConfig::SetGlobalObservables(const RooArgSet &set)
{

   if (!SetHasOnlyParameters(set, "ModelConfig::SetGlobalObservables"))
      return;

   // make global observables constant
   for (auto *arg : set) {
      arg->setAttribute("Constant", true);
   }

   fGlobalObsName = std::string(GetName()) + "_GlobalObservables";
   DefineSetInWS(fGlobalObsName.c_str(), set);
}

namespace {

std::unique_ptr<RooLinkedList>
finalizeCmdList(ModelConfig const &modelConfig, RooLinkedList const &cmdList, std::vector<RooCmdArg> &cmdArgs)
{
   auto addCmdArg = [&](RooCmdArg const &cmdArg) {
      if (cmdList.FindObject(cmdArg.GetName())) {
         std::stringstream ss;
         ss << "Illegal command argument \"" << cmdArg.GetName()
            << "\" passed to ModelConfig::createNLL(). This option is retrieved from the ModelConfig itself.";
         const std::string errorMsg = ss.str();
         oocoutE(&modelConfig, InputArguments) << errorMsg << std::endl;
         throw std::runtime_error(errorMsg);
      }
      cmdArgs.push_back(cmdArg);
   };

   if (auto args = modelConfig.GetConditionalObservables()) {
      addCmdArg(RooFit::ConditionalObservables(*args));
   }

   if (auto args = modelConfig.GetGlobalObservables()) {
      addCmdArg(RooFit::GlobalObservables(*args));
   }

   // The constraint terms of the model, and the ones that are passed in addition, e.g. a prior. They are passed also
   // if there are none: this says that the likelihood has no constraint terms, also if another ModelConfig of the
   // same pdf has some.
   RooArgSet constraints;
   if (auto modelConstraints = modelConfig.GetExternalConstraints()) {
      constraints.add(*modelConstraints);
   }
   auto finalCmdList = std::make_unique<RooLinkedList>();
   for (TObject *obj : cmdList) {
      auto *cmdArg = static_cast<RooCmdArg *>(obj);
      if (std::string{cmdArg->GetName()} != "ExternalConstraints") {
         finalCmdList->Add(cmdArg);
      } else if (RooArgSet const *extra = cmdArg->getSet(0)) {
         constraints.add(*extra, /*silent=*/true);
      }
   }
   cmdArgs.push_back(RooFit::ExternalConstraints(constraints));

   for (RooCmdArg &arg : cmdArgs) {
      finalCmdList->Add(&arg);
   }

   return finalCmdList;
}

} // namespace

/** @fn RooStats::ModelConfig::createNLL()
 *
 * Wrapper around RooAbsPdf::createNLL(), where
 * the pdf and some configuration options are retrieved from the ModelConfig.
 *
 * The options taken from the ModelConfig are:
 *
 *   * ConditionalObservables()
 *   * GlobalObservables()
 *   * ExternalConstraints()
 *
 * Constraint terms that are passed with ExternalConstraints() are added to the
 * ones of the ModelConfig, e.g. a prior pdf. Except for the other options
 * above, you can still pass all the other command arguments supported by
 * RooAbsPdf::createNLL().
 */

std::unique_ptr<RooAbsReal> ModelConfig::createNLLImpl(RooAbsData &data, const RooLinkedList &cmdList) const
{
   checkConstraintTerms();
   std::vector<RooCmdArg> cmdArgs;
   auto finalCmdList = finalizeCmdList(*this, cmdList, cmdArgs);
   return std::unique_ptr<RooAbsReal>{GetPdf()->createNLL(data, *finalCmdList)};
}

/** @fn RooStats::ModelConfig::fitTo()
 *
 * Wrapper around RooAbsPdf::fitTo(), where
 * the pdf and some configuration options are retrieved from the ModelConfig.
 *
 * See ModelConfig::createNLL() for more information.
 */
std::unique_ptr<RooFitResult> ModelConfig::fitToImpl(RooAbsData &data, const RooLinkedList &cmdList) const
{
   checkConstraintTerms();
   std::vector<RooCmdArg> cmdArgs;
   auto finalCmdList = finalizeCmdList(*this, cmdList, cmdArgs);
   return std::unique_ptr<RooFitResult>{GetPdf()->fitTo(data, *finalCmdList)};
}

} // end namespace RooStats
