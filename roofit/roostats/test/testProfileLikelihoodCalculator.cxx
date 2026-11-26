// Tests for the ProfileLikelihoodCalculator
// Authors: Jonas Rembser, CERN 2026

#include <RooStats/LikelihoodInterval.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/ProfileLikelihoodCalculator.h>

#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooHelpers.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

namespace {

// Two models of the observable x: "constrained" has the parameters mu and theta, and a constraint term of theta;
// "unconstrained" has only mu.
class PLCModels : public testing::Test {
protected:
   void SetUp() override
   {
      ws.factory("Gaussian::model(x[-10, 10], expr::mean('mu + theta', mu[0, -5, 5], theta[0, -5, 5]), 1.0)");
      ws.factory("Gaussian::model2(x, mu, 1.0)");
      ws.factory("Gaussian::constraint(gTheta[0, -5, 5], theta, 0.1)");
      ws.var("gTheta")->setConstant(true);
      addModelConfig("constrained", "model", true);
      addModelConfig("unconstrained", "model2", false);
      data = std::unique_ptr<RooDataSet>{ws.pdf("model")->generate(*ws.var("x"), 20)};
   }

   void addModelConfig(const char *name, const char *pdf, bool constrained)
   {
      RooStats::ModelConfig mc{name, &ws};
      mc.SetPdf(pdf);
      mc.SetObservables("x");
      mc.SetParametersOfInterest("mu");
      if (constrained) {
         mc.SetNuisanceParameters("theta");
         mc.SetGlobalObservables("gTheta");
         mc.SetExternalConstraints("constraint");
      }
      ws.import(mc);
   }

   RooStats::ModelConfig &modelConfig(const char *name) { return static_cast<RooStats::ModelConfig &>(*ws.obj(name)); }

   RooHelpers::LocalChangeMsgLevel changeMsgLvl{RooFit::WARNING};
   RooWorkspace ws;
   std::unique_ptr<RooDataSet> data;
};

} // namespace

// The likelihood of the calculator is the one of the ModelConfig, with its constraint terms. Without the constraint
// term, which fixes theta to about 0 within 0.1, mu would not be determined, as the data only depend on mu + theta.
TEST_F(PLCModels, LikelihoodOfTheModel)
{
   RooStats::ProfileLikelihoodCalculator plc{*data, modelConfig("constrained")};
   std::unique_ptr<RooStats::LikelihoodInterval> interval{plc.GetInterval()};
   ASSERT_NE(interval, nullptr);
   EXPECT_LT(interval->UpperLimit(*ws.var("mu")) - interval->LowerLimit(*ws.var("mu")), 1.0);
   std::unique_ptr<RooArgSet> variables{interval->GetLikelihoodRatio()->getVariables()};
   EXPECT_NE(variables->find("gTheta"), nullptr);

   // The calculator keeps no copies of the model: another model replaces everything, also the constraint terms.
   plc.SetModel(modelConfig("unconstrained"));
   interval.reset(plc.GetInterval());
   ASSERT_NE(interval, nullptr);
   variables.reset(interval->GetLikelihoodRatio()->getVariables());
   EXPECT_EQ(variables->find("gTheta"), nullptr);
   EXPECT_EQ(variables->find("theta"), nullptr);
}
