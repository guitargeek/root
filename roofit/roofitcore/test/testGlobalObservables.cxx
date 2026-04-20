// Tests for global observables
// Authors: Jonas Rembser, CERN  08/2021

#include <RooAbsPdf.h>
#include <RooCmdConfig.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHelpers.h>
#include <RooLinkedList.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRandom.h>

#include "../src/FitHelpers.h"

#include "gtest_wrapper.h"

#include <memory>
#include <functional>

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
   GlobsTest() : _evalBackend{RooFit::EvalBackend::Legacy()} {}

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

TEST_P(GlobsTest, GenerateWithGlobalObservables)
{
   using namespace RooFit;

   auto &gm = *ws().var("gm");
   auto &gs = *ws().var("gs");
   auto &x = *ws().var("x");

   // Case 1: the global observable values are taken directly from the
   // RooAbsReal instances passed to GlobalObservables().
   resetParameters();
   gm.setVal(12.5);
   gs.setVal(0.8);
   std::unique_ptr<RooDataSet> toy1{modelc().generate(x, 100, GlobalObservables(gm, gs))};

   ASSERT_NE(toy1, nullptr);
   EXPECT_EQ(toy1->numEntries(), 100);

   // The event observables should be just x.
   EXPECT_EQ(toy1->get()->size(), std::size_t(1));
   EXPECT_NE(toy1->get()->find("x"), nullptr);
   EXPECT_EQ(toy1->get()->find("gm"), nullptr);
   EXPECT_EQ(toy1->get()->find("gs"), nullptr);

   // The global observables should be stored with the current workspace
   // values.
   auto const *storedGlobs1 = toy1->getGlobalObservables();
   ASSERT_NE(storedGlobs1, nullptr);
   EXPECT_EQ(storedGlobs1->size(), std::size_t(2));
   auto const *storedGm1 = static_cast<RooRealVar const *>(storedGlobs1->find("gm"));
   auto const *storedGs1 = static_cast<RooRealVar const *>(storedGlobs1->find("gs"));
   ASSERT_NE(storedGm1, nullptr);
   ASSERT_NE(storedGs1, nullptr);
   EXPECT_DOUBLE_EQ(storedGm1->getVal(), 12.5);
   EXPECT_DOUBLE_EQ(storedGs1->getVal(), 0.8);

   // Case 2: the global observables are also passed in the set of observables
   // to generate, so their values are sampled once from the model.
   resetParameters();
   gm.setVal(12.5);
   gs.setVal(0.8);
   std::unique_ptr<RooDataSet> toy2{modelc().generate({x, gm, gs}, 100, GlobalObservables(gm, gs))};

   ASSERT_NE(toy2, nullptr);
   EXPECT_EQ(toy2->numEntries(), 100);

   // The sampled global observables should not appear as event columns.
   EXPECT_EQ(toy2->get()->size(), std::size_t(1));
   EXPECT_NE(toy2->get()->find("x"), nullptr);
   EXPECT_EQ(toy2->get()->find("gm"), nullptr);
   EXPECT_EQ(toy2->get()->find("gs"), nullptr);

   auto const *storedGlobs2 = toy2->getGlobalObservables();
   ASSERT_NE(storedGlobs2, nullptr);
   EXPECT_EQ(storedGlobs2->size(), std::size_t(2));
   auto const *storedGm2 = static_cast<RooRealVar const *>(storedGlobs2->find("gm"));
   auto const *storedGs2 = static_cast<RooRealVar const *>(storedGlobs2->find("gs"));
   ASSERT_NE(storedGm2, nullptr);
   ASSERT_NE(storedGs2, nullptr);

   // The sampled values must be within the valid range of the variables.
   EXPECT_GE(storedGm2->getVal(), gm.getMin());
   EXPECT_LE(storedGm2->getVal(), gm.getMax());
   EXPECT_GE(storedGs2->getVal(), gs.getMin());
   EXPECT_LE(storedGs2->getVal(), gs.getMax());

   // The sampled values should differ from the pre-generation values,
   // confirming that they were actually drawn from the model rather than
   // just snapshotted from the RooAbsReal instances.
   EXPECT_NE(storedGm2->getVal(), 12.5)
      << "global observable gm was not sampled when passed to the set of observables to generate";
   EXPECT_NE(storedGs2->getVal(), 0.8)
      << "global observable gs was not sampled when passed to the set of observables to generate";

   // generate() must not mutate the pdf: the workspace global observables
   // should still hold their pre-call values.
   EXPECT_DOUBLE_EQ(gm.getVal(), 12.5)
      << "generate() mutated workspace variable gm; its value should be restored after the call";
   EXPECT_DOUBLE_EQ(gs.getVal(), 0.8)
      << "generate() mutated workspace variable gs; its value should be restored after the call";

   // The stored global observables should be usable by a fit: the fit result
   // when feeding the toy back into the constrained model must be identical
   // to one where the same values are passed via a GlobalObservables command
   // argument.
   resetParameters();
   auto resFromData = doFit(modelc(), *toy2);

   resetParameters();
   gm.setVal(storedGm2->getVal());
   gs.setVal(storedGs2->getVal());
   auto resFromArgs = doFit(modelc(), *toy2, GlobalObservables(gm, gs), GlobalObservablesSource("model"));

   EXPECT_TRUE(resFromData->isIdentical(*resFromArgs))
      << "fit result using global observables stored in the generated dataset "
         "differs from the fit result using the same values passed via "
         "GlobalObservables() with source=model";
}

INSTANTIATE_TEST_SUITE_P(TestGlobalObservables, GlobsTest, testing::Values(ROOFIT_EVAL_BACKENDS),
                         [](testing::TestParamInfo<GlobsTest::ParamType> const &paramInfo) {
                            std::stringstream ss;
                            ss << "EvalBackend" << std::get<0>(paramInfo.param).name();
                            return ss.str();
                         });
