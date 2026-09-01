// Tests for global observables
// Authors: Jonas Rembser, CERN  08/2021

#include <RooAbsPdf.h>
#include <RooCmdConfig.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHelpers.h>
#include <RooLinkedList.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRandom.h>

#include <TFile.h>
#include <TSystem.h>

#include "../src/FitHelpers.h"

#include "gtest_wrapper.h"

#include <cmath>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

using RooFit::FitHelpers::minimize;

namespace {

// Helper function to check if two RooFitResults are not identical.
// We can't use RooFitResult::isIdentical() here, because it will print
// something when the comparison fails even with verbose set to false.
bool isNotIdentical(RooFitResult const &res1, RooFitResult const &res2)
{
   std::size_t n = res1.floatParsFinal().size();
   if (n != res2.floatParsFinal().size()) {
      return true;
   }
   for (std::size_t i = 0; i < n; ++i) {
      if (static_cast<RooAbsRealLValue &>(res1.floatParsFinal()[i]).getVal() !=
          static_cast<RooAbsRealLValue &>(res2.floatParsFinal()[i]).getVal())
         return true;
   }
   return false;
}

} // namespace

// Test environment to verify that if we use the feature of storing global
// observables in a RooDataSet, we can reproduce the same fit results as when
// we track the global observables separately.
class GlobsTest : public testing::TestWithParam<std::tuple<RooFit::EvalBackend>> {
public:
   GlobsTest() : _evalBackend{RooFit::EvalBackend(RooFit::EvalBackend::Value::Legacy)} {}

   void SetUp() override
   {
      RooRandom::randomGenerator()->SetSeed(1337ul);

      // silence log output
      _changeMsgLvl = std::make_unique<RooHelpers::LocalChangeMsgLevel>(RooFit::WARNING);

      _evalBackend = std::get<0>(GetParam());

      // We use the global observable also in the model for the event
      // observables. It's unusual, but let's better do this to also cover the
      // corner case where the global observable is not only part of the
      // constraint term.
      _ws.factory("Product::sigma({s[2.0, 0.1, 10.0], gs[1.0, 0.1, 10.0]})");

      _ws.factory("Gaussian::model(x[0.0, 0.0, 20.0], m[10.0, 0.0, 20.0], sigma)");

      // the constraint pdfs, they are RooPoisson so we can't have tests that accidentally
      // pass because of the symmetry of normalizing over x or mu
      _ws.factory("Poisson::mconstraint(gm[11.0, 0.0, 20.0], m)");
      _ws.factory("Poisson::sconstraint(gs, s)");

      // global observables, always constant in fits
      RooRealVar &gm = *_ws.var("gm");
      RooRealVar &gs = *_ws.var("gs");
      gm.setConstant(true);
      gs.setConstant(true);

      // the model multiplied with the constraint term
      _ws.factory("ProdPdf::modelc({model, mconstraint, sconstraint})");

      // generate small dataset for use in fitting below, also cloned versions
      // with one or two global observables attached
      _data = std::unique_ptr<RooDataSet>{_ws.pdf("model")->generate(*_ws.var("x"), 50)};

      _dataWithMeanSigmaGlobs.reset(static_cast<RooDataSet *>(_data->Clone()));
      _dataWithMeanSigmaGlobs->SetName((std::string(_data->GetName()) + "_gm_gs").c_str());
      _dataWithMeanSigmaGlobs->setGlobalObservables({gm, gs});

      _dataWithMeanGlob.reset(static_cast<RooDataSet *>(_data->Clone((std::string(_data->GetName()) + "_gm").c_str())));
      _dataWithMeanGlob->setGlobalObservables(gm);
   }

   // reset the parameter values to initial values before fits
   void resetParameters()
   {
      std::vector<std::string> names{"x", "m", "s", "gm", "gs"};
      std::vector<double> values{0.0, 10.0, 2.0, 11.0, 1.0};
      for (std::size_t i = 0; i < names.size(); ++i) {
         auto *var = _ws.var(names[i]);
         var->setVal(values[i]);
         var->setError(0.0);
      }
   }

   RooFit::EvalBackend const &evalBackend() { return _evalBackend; }
   RooWorkspace &ws() { return _ws; }
   RooDataSet &data() { return *_data; }
   RooDataSet &dataWithMeanSigmaGlobs() { return *_dataWithMeanSigmaGlobs; }
   RooDataSet &dataWithMeanGlob() { return *_dataWithMeanGlob; }
   RooAbsPdf &model() { return *ws().pdf("model"); }
   RooAbsPdf &modelc() { return *ws().pdf("modelc"); }

   std::unique_ptr<RooFitResult> doFit(RooAbsPdf &model, RooAbsData &data, RooCmdArg const &arg1 = {},
                                       RooCmdArg const &arg2 = {}, RooCmdArg const &arg3 = {},
                                       RooCmdArg const &arg4 = {})
   {
      using namespace RooFit;
      return std::unique_ptr<RooFitResult>{
         model.fitTo(data, Save(), Verbose(false), PrintLevel(-1), _evalBackend, arg1, arg2, arg3, arg4)};
   }

   void TearDown() override
   {
      _data.reset();
      _dataWithMeanSigmaGlobs.reset();
      _data.reset();
      _changeMsgLvl.reset();
   }

private:
   RooFit::EvalBackend _evalBackend;
   RooWorkspace _ws;
   std::unique_ptr<RooDataSet> _data;
   std::unique_ptr<RooDataSet> _dataWithMeanSigmaGlobs;
   std::unique_ptr<RooDataSet> _dataWithMeanGlob;
   std::unique_ptr<RooHelpers::LocalChangeMsgLevel> _changeMsgLvl;
};

TEST_P(GlobsTest, NoConstraints)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");

   // fit with no constraints
   resetParameters();
   auto res1 = doFit(model(), data());
   resetParameters();
   // vary global observable to verify true value is picked up from the dataset
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   double gmVaryVal = gm.getVal();
   double gsVaryVal = gs.getVal();
   auto res2 = doFit(model(), dataWithMeanSigmaGlobs());
   EXPECT_TRUE(res1->isIdentical(*res2)) << "fitting an unconstrained model "
                                            "gave a different result when unrelated global observables were stored in "
                                            "the dataset";

   // verify that taking the global observable values from data has not changed
   // the values in the model
   {
      const auto message = "taking the global observable values from data has changed the values in the model";
      EXPECT_EQ(gmVaryVal, gm.getVal()) << message;
      EXPECT_EQ(gsVaryVal, gs.getVal()) << message;
   }
}

TEST_P(GlobsTest, InternalConstraints)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");

   // constrained fit with RooProdPdf
   resetParameters();
   auto res1 = doFit(modelc(), data(), GlobalObservables(gm, gs));
   resetParameters();
   // vary global observable to verify true value is picked up from the dataset
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   double gmVaryVal = gm.getVal();
   double gsVaryVal = gs.getVal();
   auto res2 = doFit(modelc(), dataWithMeanSigmaGlobs());
   EXPECT_TRUE(res1->isIdentical(*res2)) << "fitting an model with internal "
                                            "constraints in a RooPrdPdf gave a different result when global "
                                            "observables were stored in the dataset";

   // verify that taking the global observable values from data has not changed
   // the values in the model
   {
      const auto message = "taking the global observable values from data has changed the values in the model";
      EXPECT_EQ(gmVaryVal, gm.getVal()) << message;
      EXPECT_EQ(gsVaryVal, gs.getVal()) << message;
   }
}

TEST_P(GlobsTest, ExternalConstraints)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");
   auto &mconstraint = *ws().pdf("mconstraint");
   auto &sconstraint = *ws().pdf("sconstraint");

   // constrained fit with external constraints
   resetParameters();
   auto res1 = doFit(model(), data(), ExternalConstraints({mconstraint, sconstraint}), GlobalObservables(gm, gs));
   resetParameters();
   // vary global observable to verify true value is picked up from the dataset
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   double gmVaryVal = gm.getVal();
   double gsVaryVal = gs.getVal();
   auto res2 = doFit(model(), dataWithMeanSigmaGlobs(), ExternalConstraints({mconstraint, sconstraint}));
   EXPECT_TRUE(res1->isIdentical(*res2))
      << "fitting an model with external "
         "constraints passed via ExternalConstraints() gave a different result when global "
         "observables were stored in the dataset";

   // verify that taking the global observable values from data has not changed
   // the values in the model
   {
      const auto message = "taking the global observable values from data has changed the values in the model";
      EXPECT_EQ(gmVaryVal, gm.getVal()) << message;
      EXPECT_EQ(gsVaryVal, gs.getVal()) << message;
   }
}

TEST_P(GlobsTest, SubsetOfConstraintsFromData)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");

   // check if only a subset of constraints it taken from data
   resetParameters();
   auto res1 = doFit(modelc(), data(), GlobalObservables(gm, gs));
   resetParameters();
   // Vary global observable to verify true value is picked up from the dataset.
   // This time we only get gm from the dataset, so we don't vary gs for now.
   gm.setVal(gm.getVal() + 0.5);
   double gmVaryVal = gm.getVal();
   double gsVaryVal = gs.getVal();
   // if we take the global observables from the model, they have to be constant:
   auto res2 = doFit(modelc(), dataWithMeanGlob(), GlobalObservables(gm, gs));
   EXPECT_TRUE(res1->isIdentical(*res2)) << "fitting a constrained model "
                                            "to a dataset that only stores a subset of the defined global observables "
                                            "gave the wrong result";

   // verify that taking the global observable values from data has not changed
   // the values in the model
   {
      const auto message = "taking the global observable values from data has changed the values in the model";
      EXPECT_EQ(gmVaryVal, gm.getVal()) << message;
      EXPECT_EQ(gsVaryVal, gs.getVal()) << message;
   }

   resetParameters();
   // Now that we also vary gs, the fit results should not be identical.
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   gmVaryVal = gm.getVal();
   gsVaryVal = gs.getVal();
   auto res3 = doFit(modelc(), dataWithMeanGlob(), GlobalObservables(gm, gs));
   EXPECT_TRUE(isNotIdentical(*res1, *res3))
      << "fitting a constrained model "
         "to a dataset that only stores a subset of the defined global observables "
         "gave the wrong result";

   // verify that taking the global observable values from data has not changed
   // the values in the model
   {
      const auto message = "taking the global observable values from data has changed the values in the model";
      EXPECT_EQ(gmVaryVal, gm.getVal()) << message;
      EXPECT_EQ(gsVaryVal, gs.getVal()) << message;
   }
}

namespace {

RooCmdConfig minimizerCfg()
{

   RooCmdConfig pc("minimizerCfg");

   RooFit::FitHelpers::defineMinimizationOptions(pc);

   std::vector<RooCmdArg> cmdArgs{RooFit::Save(), RooFit::PrintLevel(-1)};

   RooLinkedList cmdList;
   for (auto &arg : cmdArgs) {
      cmdList.Add(&arg);
   }

   pc.process(cmdList);

   return pc;
}

} // namespace

TEST_P(GlobsTest, ResetDataToWrongData)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");
   auto &model = modelc();

   // constrained fit with RooProdPdf
   resetParameters();
   auto res1 = doFit(model, data(), GlobalObservables(gm, gs));

   resetParameters();
   // vary global observable to deliberately store "wrong" values in a cloned dataset
   std::unique_ptr<RooDataSet> wrongData{static_cast<RooDataSet *>(data().Clone())};
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   wrongData->setGlobalObservables({gm, gs});

   // check that the fit works when using the dataset with the correct values
   std::unique_ptr<RooAbsReal> nll{model.createNLL(dataWithMeanSigmaGlobs(), EvalBackend(evalBackend()))};
   auto res2 = minimize(model, *nll, dataWithMeanSigmaGlobs(), minimizerCfg());
   EXPECT_TRUE(res1->isIdentical(*res2)) << "fitting an model with internal "
                                            "constraints in a RooPrdPdf gave a different result when global "
                                            "observables were stored in the dataset";

   nll->setData(*wrongData);
   resetParameters();
   auto res3 = minimize(model, *nll, *wrongData, minimizerCfg());

   // If resetting the dataset used for the nll worked correctly also for
   // global observables, the fit will now give the wrong result.
   EXPECT_TRUE(isNotIdentical(*res1, *res3))
      << "resetting the dataset "
         "underlying a RooNLLVar didn't change the global observable value, but it "
         "should have";
}

TEST_P(GlobsTest, ResetDataToCorrectData)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");
   auto &model = modelc();

   // constrained fit with RooProdPdf
   resetParameters();
   auto res1 = doFit(model, data(), GlobalObservables(gm, gs));

   resetParameters();
   // vary global observable to deliberately store "wrong" values in a cloned dataset
   std::unique_ptr<RooDataSet> wrongData{static_cast<RooDataSet *>(data().Clone())};
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   wrongData->setGlobalObservables({gm, gs});
   resetParameters();

   // check that the fit doesn't work when using the dataset with the wrong values
   std::unique_ptr<RooAbsReal> nll{model.createNLL(*wrongData, EvalBackend(evalBackend()))};
   auto res2 = minimize(model, *nll, *wrongData, minimizerCfg());
   EXPECT_TRUE(isNotIdentical(*res1, *res2)) << "fitting an model with internal "
                                                "constraints in a RooPrdPdf ignored the global "
                                                "observables stored in the dataset";

   nll->setData(dataWithMeanSigmaGlobs());
   resetParameters();
   auto res3 = minimize(model, *nll, dataWithMeanSigmaGlobs(), minimizerCfg());
   EXPECT_TRUE(res1->isIdentical(*res3)) << "resetting the dataset "
                                            "underlying a RooNLLVar didn't change the global observable value, but it "
                                            "should have";
}

TEST_P(GlobsTest, GlobalObservablesSourceFromModel)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");

   // constrained fit with RooProdPdf
   resetParameters();
   auto res1 = doFit(modelc(), data(), GlobalObservables(gm, gs));
   resetParameters();
   // vary global observable to verify true value is picked up from the dataset
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);

   // verify that fit results are identical when global observable values are
   // taken from data
   auto res2 = doFit(modelc(), dataWithMeanSigmaGlobs());
   EXPECT_TRUE(res1->isIdentical(*res2));

   auto res3 = doFit(modelc(), dataWithMeanSigmaGlobs(), GlobalObservablesSource("model"), GlobalObservables(gm, gs));

   // If the global observable values are indeed taken from the model and not
   // from data, the comparison will fail now because we have changed the
   // global observable values of the model after the first fit.
   EXPECT_TRUE(isNotIdentical(*res2, *res3));
}

TEST_P(GlobsTest, ResetDataButSourceFromModel)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");
   auto &model = modelc();

   // constrained fit with RooProdPdf
   resetParameters();
   auto res1 = doFit(model, data(), GlobalObservables(gm, gs));

   resetParameters();
   // vary global observable to deliberately store "wrong" values in a cloned dataset
   std::unique_ptr<RooDataSet> wrongData{static_cast<RooDataSet *>(data().Clone())};
   gm.setVal(gm.getVal() + 0.5);
   gs.setVal(gs.getVal() + 2.5);
   wrongData->setGlobalObservables({gm, gs});

   resetParameters();

   // check that the fit works when using the dataset with the correct values
   std::unique_ptr<RooAbsReal> nll{model.createNLL(dataWithMeanSigmaGlobs(), GlobalObservablesSource("model"),
                                                   GlobalObservables(gm, gs), EvalBackend(evalBackend()))};
   auto res2 = minimize(model, *nll, dataWithMeanSigmaGlobs(), minimizerCfg());
   EXPECT_TRUE(res1->isIdentical(*res2));

   nll->setData(*wrongData);
   resetParameters();
   auto res3 = minimize(model, *nll, *wrongData, minimizerCfg());

   // this time it should still be identical because even though we reset to
   // the wrong data, we set the global observables source to "model"
   EXPECT_TRUE(res1->isIdentical(*res3));
}

INSTANTIATE_TEST_SUITE_P(TestGlobalObservables, GlobsTest, testing::Values(ROOFIT_EVAL_BACKENDS),
                         [](testing::TestParamInfo<GlobsTest::ParamType> const &paramInfo) {
                            std::stringstream ss;
                            ss << "EvalBackend" << std::get<0>(paramInfo.param).name();
                            return ss.str();
                         });

////////////////////////////////////////////////////////////////////////////////
/// Tests for the global observables support in toy dataset generation
/// (GitHub issue #10634).

namespace {

/// Model with one event observable `x` and one global observable `gmu` that
/// constrains the nuisance parameter `mu` with a Gaussian of width `sigmaG`.
std::unique_ptr<RooWorkspace> makeConstrainedModel()
{
   auto ws = std::make_unique<RooWorkspace>("ws");
   ws->factory("Gaussian::model(x[-10, 10], mu[0.0, -10, 10], sigma[2.0, 0.1, 10.0])");
   ws->factory("Gaussian::constraint(gmu[0.0, -10, 10], mu, sigmaG[1.5, 0.01, 10.0])");
   ws->factory("ProdPdf::modelc({model, constraint})");
   ws->var("gmu")->setConstant(true);
   return ws;
}

double globValue(RooAbsData const &data, const char *name)
{
   RooArgSet const *globs = data.getGlobalObservables();
   if (!globs || !globs->find(name)) {
      throw std::runtime_error(std::string("no global observable ") + name + " in dataset!");
   }
   return static_cast<RooRealVar const *>(globs->find(name))->getVal();
}

/// Sample mean and (unbiased) sample standard deviation of a set of toy values.
struct Moments {
   double mean = 0.0;
   double stddev = 0.0;
};

Moments moments(std::vector<double> const &values)
{
   const double n = values.size();
   double sum = 0.0;
   double sum2 = 0.0;
   for (double v : values) {
      sum += v;
      sum2 += v * v;
   }
   Moments out;
   out.mean = sum / n;
   out.stddev = std::sqrt((sum2 - n * out.mean * out.mean) / (n - 1.));
   return out;
}

/// Generate `nToys` datasets from `spec` and return the sampled values of the
/// global observable `name`.
std::vector<double> sampleGlobs(RooAbsPdf &pdf, RooAbsPdf::GenSpec &spec, const char *name, int nToys)
{
   std::vector<double> out;
   out.reserve(nToys);
   for (int i = 0; i < nToys; ++i) {
      std::unique_ptr<RooDataSet> data{pdf.generate(spec)};
      out.push_back(globValue(*data, name));
   }
   return out;
}

} // namespace

/// Mode 1: the global observable is not in the set of variables to generate,
/// so its current value is taken from the model and stored in the dataset.
TEST(GlobalObservablesGeneration, TakeValueFromModel)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   gmu.setVal(1.234);

   std::unique_ptr<RooDataSet> data{modelc.generate(x, 100, RooFit::GlobalObservables(gmu))};

   ASSERT_TRUE(data != nullptr);
   EXPECT_EQ(data->numEntries(), 100);

   // The global observable is not a column of the dataset
   EXPECT_EQ(data->get()->find("gmu"), nullptr);

   RooArgSet const *globs = data->getGlobalObservables();
   ASSERT_TRUE(globs != nullptr);
   EXPECT_EQ(globs->size(), 1u);
   EXPECT_DOUBLE_EQ(globValue(*data, "gmu"), 1.234);
   // Global observables attached to a dataset are always constant
   EXPECT_TRUE(static_cast<RooRealVar const *>(globs->find("gmu"))->isConstant());

   // The dataset owns a snapshot: changing the model doesn't change the dataset
   gmu.setVal(5.0);
   EXPECT_DOUBLE_EQ(globValue(*data, "gmu"), 1.234);
}

/// Mode 2: the global observable is also in the set of variables to generate,
/// so its value is sampled from the constraint term in the model. Check that
/// the sampled values follow the constraint pdf.
TEST(GlobalObservablesGeneration, SampleFromConstraint)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &mu = *ws->var("mu");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   // The global observable is sampled from constraint(gmu | mu, sigmaG) with
   // the *current* value of the nuisance parameter mu.
   mu.setVal(1.0);
   const double muRef = mu.getVal();
   const double sigmaGRef = ws->var("sigmaG")->getVal();

   // Set the global observable to a value far away from the constraint pdf, so
   // we can verify that it is really sampled and that the model is untouched.
   gmu.setVal(-9.0);

   // Using the prepareMultiGen()/generate() combination also covers the
   // GenSpec code path, and it is much faster for many toys.
   std::unique_ptr<RooAbsPdf::GenSpec> spec{
      modelc.prepareMultiGen({x, gmu}, RooFit::NumEvents(10), RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(spec != nullptr);

   const int nToys = 5000;

   // Structural checks on the first toy
   {
      std::unique_ptr<RooDataSet> data{modelc.generate(*spec)};
      ASSERT_TRUE(data != nullptr);
      EXPECT_EQ(data->numEntries(), 10);
      // Sampled global observables are not columns of the dataset: there is one
      // auxiliary measurement per toy experiment, not one per event.
      EXPECT_EQ(data->get()->find("gmu"), nullptr);
   }

   std::vector<double> values = sampleGlobs(modelc, *spec, "gmu", nToys);
   const Moments mom = moments(values);

   // Four times the statistical uncertainty on the mean and on the standard
   // deviation. This is a tight check: it fails if the values are not sampled
   // at all (they would all be -9.0), if they are sampled from the wrong pdf
   // (e.g. from model(x | mu, sigma), which has a different width), or if they
   // are sampled uniformly over the range of gmu.
   EXPECT_NEAR(mom.mean, muRef, 4.0 * sigmaGRef / std::sqrt(double(nToys)));
   EXPECT_NEAR(mom.stddev, sigmaGRef, 4.0 * sigmaGRef / std::sqrt(2.0 * nToys));

   // Each toy gets its own auxiliary measurement, so the values must not repeat
   EXPECT_EQ(std::set<double>(values.begin(), values.end()).size(), values.size());

   // The mean of the sampled values follows the nuisance parameter: sampling
   // again with a shifted mu must shift the distribution by exactly that much.
   // This is what distinguishes sampling from the constraint term from sampling
   // from any other distribution that happens to be centered around muRef.
   const double muShifted = -2.75;
   mu.setVal(muShifted);
   const Moments momShifted = moments(sampleGlobs(modelc, *spec, "gmu", nToys));
   EXPECT_NEAR(momShifted.mean, muShifted, 4.0 * sigmaGRef / std::sqrt(double(nToys)));
   EXPECT_NEAR(momShifted.stddev, sigmaGRef, 4.0 * sigmaGRef / std::sqrt(2.0 * nToys));
   mu.setVal(muRef);

   // The width of the sampled values follows the width of the constraint pdf
   const double sigmaGShifted = 0.4;
   ws->var("sigmaG")->setVal(sigmaGShifted);
   const Moments momNarrow = moments(sampleGlobs(modelc, *spec, "gmu", nToys));
   EXPECT_NEAR(momNarrow.mean, muRef, 4.0 * sigmaGShifted / std::sqrt(double(nToys)));
   EXPECT_NEAR(momNarrow.stddev, sigmaGShifted, 4.0 * sigmaGShifted / std::sqrt(2.0 * nToys));
   ws->var("sigmaG")->setVal(sigmaGRef);

   // The state of the model is not changed by the sampling
   EXPECT_DOUBLE_EQ(gmu.getVal(), -9.0);

   // The same works with the plain generate() interface
   std::unique_ptr<RooDataSet> data{modelc.generate({x, gmu}, 10, RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(data != nullptr);
   EXPECT_EQ(data->get()->find("gmu"), nullptr);
   EXPECT_NE(globValue(*data, "gmu"), -9.0);
   EXPECT_DOUBLE_EQ(gmu.getVal(), -9.0);
}

/// The generated dataset is consumed correctly by createNLL(), which picks up
/// the global observables from the dataset by default.
TEST(GlobalObservablesGeneration, RoundTripToNll)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   gmu.setVal(1.0);
   std::unique_ptr<RooDataSet> data{modelc.generate(x, 100, RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(data != nullptr);

   const double nllRef = std::unique_ptr<RooAbsReal>{modelc.createNLL(*data)}->getVal();

   // Changing the global observable in the model doesn't change the NLL,
   // because the value is taken from the dataset by default
   gmu.setVal(3.0);
   EXPECT_DOUBLE_EQ(std::unique_ptr<RooAbsReal>{modelc.createNLL(*data)}->getVal(), nllRef);

   // ... unless the values are explicitly requested to come from the model
   const double nllFromModel = std::unique_ptr<RooAbsReal>{modelc.createNLL(*data, RooFit::GlobalObservables(gmu),
                                                                            RooFit::GlobalObservablesSource("model"))}
                                  ->getVal();
   EXPECT_NE(nllFromModel, nllRef);

   // The difference is exactly the change in the Gaussian constraint term
   const double sigmaG = ws->var("sigmaG")->getVal();
   const double mu = ws->var("mu")->getVal();
   auto constraintNll = [&](double g) { return 0.5 * (g - mu) * (g - mu) / (sigmaG * sigmaG); };
   EXPECT_NEAR(nllFromModel - nllRef, constraintNll(3.0) - constraintNll(1.0), 1e-9);
}

/// Error cases.
TEST(GlobalObservablesGeneration, ErrorCases)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::FATAL};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &model = *ws->pdf("model");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   // The model doesn't depend on the global observable, so there is no
   // constraint term from which it could be sampled
   EXPECT_THROW(std::unique_ptr<RooDataSet>{model.generate({x, gmu}, 10, RooFit::GlobalObservables(gmu))},
                std::invalid_argument);
   EXPECT_THROW(std::unique_ptr<RooAbsPdf::GenSpec>{model.prepareMultiGen({x, gmu}, RooFit::NumEvents(10),
                                                                          RooFit::GlobalObservables(gmu))},
                std::invalid_argument);

   // Nothing left to generate: the dataset would have no columns
   EXPECT_THROW(std::unique_ptr<RooDataSet>{modelc.generate(gmu, 10, RooFit::GlobalObservables(gmu))},
                std::invalid_argument);

   // Storing the value of a global observable that is not part of the model is
   // allowed: no sampling is involved
   RooRealVar gOther{"gOther", "gOther", 42.0};
   std::unique_ptr<RooDataSet> data{model.generate(x, 10, RooFit::GlobalObservables(gOther))};
   ASSERT_TRUE(data != nullptr);
   EXPECT_DOUBLE_EQ(globValue(*data, "gOther"), 42.0);
}

/// The global observables are also attached to binned datasets generated with
/// RooAbsPdf::generateBinned(), with the same two modes.
TEST(GlobalObservablesGeneration, GenerateBinned)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &mu = *ws->var("mu");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   // Mode 1: value taken from the model
   gmu.setVal(1.234);
   std::unique_ptr<RooDataHist> hist{modelc.generateBinned(x, 500, RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(hist != nullptr);
   EXPECT_EQ(hist->get()->find("gmu"), nullptr);
   EXPECT_DOUBLE_EQ(globValue(*hist, "gmu"), 1.234);

   // An Asimov dataset gets the global observables at their nominal values too
   std::unique_ptr<RooDataHist> asimov{
      modelc.generateBinned(x, 500, RooFit::ExpectedData(), RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(asimov != nullptr);
   EXPECT_DOUBLE_EQ(globValue(*asimov, "gmu"), 1.234);

   // Mode 2: sampled from the constraint term
   mu.setVal(-1.25);
   const double sigmaGRef = ws->var("sigmaG")->getVal();
   const int nToys = 2000;
   std::vector<double> values;
   values.reserve(nToys);
   for (int i = 0; i < nToys; ++i) {
      std::unique_ptr<RooDataHist> data{modelc.generateBinned({x, gmu}, 20, RooFit::GlobalObservables(gmu))};
      ASSERT_TRUE(data != nullptr);
      // The sampled global observable is not an axis of the histogram
      EXPECT_EQ(data->get()->find("gmu"), nullptr);
      values.push_back(globValue(*data, "gmu"));
   }
   const Moments mom = moments(values);
   EXPECT_NEAR(mom.mean, mu.getVal(), 4.0 * sigmaGRef / std::sqrt(double(nToys)));
   EXPECT_NEAR(mom.stddev, sigmaGRef, 4.0 * sigmaGRef / std::sqrt(2.0 * nToys));
   // The state of the model is not changed by the sampling
   EXPECT_DOUBLE_EQ(gmu.getVal(), 1.234);

   // Same error handling as in generate()
   EXPECT_THROW(
      std::unique_ptr<RooDataHist>{ws->pdf("model")->generateBinned({x, gmu}, 10, RooFit::GlobalObservables(gmu))},
      std::invalid_argument);
}

/// Sampling global observables also works for RooSimultaneous models, where the
/// global observables are not associated to any specific channel.
TEST(GlobalObservablesGeneration, Simultaneous)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   RooWorkspace ws{"ws"};
   ws.factory("Gaussian::pdfA(x[-10, 10], muA[1.1, -10, 10], sA[2.0])");
   ws.factory("Gaussian::pdfB(x, muB[-0.4, -10, 10], sB[1.5])");
   ws.factory("Gaussian::cA(gA[0.0, -10, 10], muA, sgA[1.0])");
   ws.factory("Gaussian::cB(gB[0.0, -10, 10], muB, sgB[0.4])");
   ws.factory("SUM::epdfA(nA[40, 0, 10000] * pdfA)");
   ws.factory("SUM::epdfB(nB[60, 0, 10000] * pdfB)");
   ws.factory("ProdPdf::mA({epdfA, cA})");
   ws.factory("ProdPdf::mB({epdfB, cB})");
   ws.factory("SIMUL::sim(cat[A=0, B=1], A=mA, B=mB)");
   ws.var("gA")->setConstant(true);
   ws.var("gB")->setConstant(true);

   RooAbsPdf &sim = *ws.pdf("sim");
   RooArgSet globs{*ws.var("gA"), *ws.var("gB")};
   RooArgSet whatVars{*ws.var("x"), *ws.cat("cat"), *ws.var("gA"), *ws.var("gB")};

   std::unique_ptr<RooAbsPdf::GenSpec> spec{
      sim.prepareMultiGen(whatVars, RooFit::Extended(), RooFit::GlobalObservables(globs))};
   ASSERT_TRUE(spec != nullptr);

   const int nToys = 2000;
   std::vector<double> valuesA;
   std::vector<double> valuesB;
   valuesA.reserve(nToys);
   valuesB.reserve(nToys);
   for (int i = 0; i < nToys; ++i) {
      std::unique_ptr<RooDataSet> data{sim.generate(*spec)};
      ASSERT_TRUE(data != nullptr);
      EXPECT_EQ(data->get()->find("gA"), nullptr);
      EXPECT_EQ(data->get()->find("gB"), nullptr);
      valuesA.push_back(globValue(*data, "gA"));
      valuesB.push_back(globValue(*data, "gB"));
   }

   // Each global observable follows the constraint of its own channel
   const Moments momA = moments(valuesA);
   const Moments momB = moments(valuesB);
   const double sgA = ws.var("sgA")->getVal();
   const double sgB = ws.var("sgB")->getVal();
   EXPECT_NEAR(momA.mean, ws.var("muA")->getVal(), 4.0 * sgA / std::sqrt(double(nToys)));
   EXPECT_NEAR(momA.stddev, sgA, 4.0 * sgA / std::sqrt(2.0 * nToys));
   EXPECT_NEAR(momB.mean, ws.var("muB")->getVal(), 4.0 * sgB / std::sqrt(double(nToys)));
   EXPECT_NEAR(momB.stddev, sgB, 4.0 * sgB / std::sqrt(2.0 * nToys));

   // The values are attached to the dataset and picked up by createNLL()
   std::unique_ptr<RooDataSet> data{sim.generate(*spec)};
   const double nllRef = std::unique_ptr<RooAbsReal>{sim.createNLL(*data)}->getVal();
   ws.var("gA")->setVal(5.0);
   ws.var("gB")->setVal(-5.0);
   EXPECT_DOUBLE_EQ(std::unique_ptr<RooAbsReal>{sim.createNLL(*data)}->getVal(), nllRef);
}

/// The global observables attached to the generated dataset survive the round
/// trip through a file, so that toys can be generated and fitted in separate
/// jobs.
TEST(GlobalObservablesGeneration, FilePersistence)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   gmu.setVal(0.77);
   std::unique_ptr<RooDataSet> data{modelc.generate(x, 100, RooFit::GlobalObservables(gmu))};
   ASSERT_TRUE(data != nullptr);
   const double nllRef = std::unique_ptr<RooAbsReal>{modelc.createNLL(*data)}->getVal();

   const char *fileName = "testGlobalObservablesGeneration.root";
   {
      TFile file{fileName, "RECREATE"};
      file.WriteObject(data.get(), "data");
   }

   {
      TFile file{fileName};
      auto *dataFromFile = file.Get<RooDataSet>("data");
      ASSERT_TRUE(dataFromFile != nullptr);
      ASSERT_TRUE(dataFromFile->getGlobalObservables() != nullptr);
      EXPECT_DOUBLE_EQ(globValue(*dataFromFile, "gmu"), 0.77);
      // The reconstructed dataset gives the same likelihood
      gmu.setVal(-3.0);
      EXPECT_DOUBLE_EQ(std::unique_ptr<RooAbsReal>{modelc.createNLL(*dataFromFile)}->getVal(), nllRef);
   }

   gSystem->Unlink(fileName);
}

/// The stored global observables actually steer a subsequent fit with default
/// arguments.
TEST(GlobalObservablesGeneration, RoundTripToFit)
{
   RooRandom::randomGenerator()->SetSeed(1337ul);
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooRealVar &mu = *ws->var("mu");
   RooRealVar &gmu = *ws->var("gmu");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   // With a constraint that is tight compared to the 100 generated events, the
   // fitted nuisance parameter is dragged towards the auxiliary measurement.
   ws->var("sigmaG")->setVal(0.3);
   ws->var("sigmaG")->setConstant(true);
   ws->var("sigma")->setConstant(true);

   auto fitMu = [&](double trueMu) {
      mu.setVal(trueMu);
      gmu.setVal(trueMu);
      std::unique_ptr<RooDataSet> data{modelc.generate(x, 100, RooFit::GlobalObservables(gmu))};
      // The value in the model is irrelevant: what counts is what is in the
      // data. If the constraint used the model value, the fit would end up
      // more than one unit away from trueMu.
      gmu.setVal(-8.0);
      mu.setVal(0.0);
      std::unique_ptr<RooFitResult> res{
         modelc.fitTo(*data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(false))};
      EXPECT_EQ(res->status(), 0);
      return static_cast<RooRealVar const *>(res->floatParsFinal().find("mu"))->getVal();
   };

   EXPECT_NEAR(fitMu(1.5), 1.5, 0.6);
   EXPECT_NEAR(fitMu(-2.0), -2.0, 0.6);
}

/// An empty GlobalObservables() set has no effect at all. In particular, no
/// empty snapshot must be attached to the dataset, because RooAbsPdf::fitTo()
/// would then normalize the constraint terms with respect to nothing.
TEST(GlobalObservablesGeneration, EmptySet)
{
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   auto ws = makeConstrainedModel();
   RooRealVar &x = *ws->var("x");
   RooAbsPdf &modelc = *ws->pdf("modelc");

   RooRandom::randomGenerator()->SetSeed(1337ul);
   std::unique_ptr<RooDataSet> ref{modelc.generate(x, 100)};
   RooRandom::randomGenerator()->SetSeed(1337ul);
   std::unique_ptr<RooDataSet> data{modelc.generate(x, 100, RooFit::GlobalObservables(RooArgSet{}))};

   ASSERT_TRUE(data != nullptr);
   EXPECT_EQ(data->getGlobalObservables(), nullptr);
   EXPECT_DOUBLE_EQ(std::unique_ptr<RooAbsReal>{modelc.createNLL(*data)}->getVal(),
                    std::unique_ptr<RooAbsReal>{modelc.createNLL(*ref)}->getVal());

   std::unique_ptr<RooAbsPdf::GenSpec> spec{
      modelc.prepareMultiGen({x}, RooFit::NumEvents(10), RooFit::GlobalObservables(RooArgSet{}))};
   ASSERT_TRUE(spec != nullptr);
   std::unique_ptr<RooDataSet> fromSpec{modelc.generate(*spec)};
   EXPECT_EQ(fromSpec->getGlobalObservables(), nullptr);
}
