// Tests for the helpers for the constraint terms of a model, see RooFit/Constraints.h
// Authors: Jonas Rembser, CERN 2026

#include <RooFit/Constraints.h>

#include <RooFit/ModelConfig.h>

#include <RooAbsPdf.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooHelpers.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

#include "gtest_wrapper.h"

#include <cmath>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

namespace {

std::set<std::string> names(RooAbsCollection const &coll)
{
   std::set<std::string> out;
   for (RooAbsArg const *arg : coll) {
      out.insert(arg->GetName());
   }
   return out;
}

using Names = std::set<std::string>;

// A model with two channels over the observable x: A has the parameters mu and alpha, B has mu and beta. The
// parameters alpha, beta and mu have constraint terms, and so has gamma, which is no parameter of the model. The
// pdf sigY is a second measurement of mu, in the observable y.
class Constraints : public testing::Test {
protected:
   void SetUp() override
   {
      for (const char *expression :
           {"Gaussian::sigA(x[-10, 10], expr::meanA('mu + alpha', mu[0, -5, 5], alpha[0, -5, 5]), 1.0)",
            "Gaussian::sigB(x, expr::meanB('mu + beta', mu, beta[0, -5, 5]), 2.0)",
            "Gaussian::sigY(y[-10, 10], mu, 1.5)", "SIMUL::model(channel[A, B], A=sigA, B=sigB)",
            "Gaussian::cAlpha(gAlpha[0, -5, 5], alpha, 1.0)", "Gaussian::cBeta(gBeta[0, -5, 5], beta, 1.0)",
            "Gaussian::cMu(gMu[0, -5, 5], mu, 3.0)", "Gaussian::cGamma(gGamma[0, -5, 5], gamma[0, -5, 5], 1.0)"}) {
         ASSERT_NE(ws.factory(expression), nullptr) << expression;
      }
      for (const char *glob : {"gAlpha", "gBeta", "gMu", "gGamma"}) {
         ws.var(glob)->setConstant(true);
      }
   }

   RooAbsPdf &pdf(const char *name)
   {
      RooAbsPdf *out = ws.pdf(name);
      if (!out) {
         throw std::runtime_error(std::string("no pdf ") + name);
      }
      return *out;
   }

   // A dataset for the channels A and B with a few values of x.
   std::unique_ptr<RooDataSet> makeData()
   {
      RooRealVar &x = *ws.var("x");
      RooCategory &channel = *ws.cat("channel");
      auto data = std::make_unique<RooDataSet>("data", "data", RooArgSet{x, channel});
      for (auto const &[label, value] :
           std::vector<std::pair<const char *, double>>{{"A", -0.4}, {"A", 0.3}, {"A", 1.1}, {"B", -1.5}, {"B", 2.0}}) {
         channel.setLabel(label);
         x.setVal(value);
         data->add({x, channel});
      }
      return data;
   }

   void addModelConfig(const char *name, const char *pdfName, RooArgSet const *constraints = nullptr)
   {
      RooFit::ModelConfig mc{name, &ws};
      mc.SetPdf(pdfName);
      mc.SetObservables("x");
      mc.SetParametersOfInterest("mu");
      if (constraints) {
         mc.SetExternalConstraints(*constraints);
      }
      ws.import(mc);
   }

   RooHelpers::LocalChangeMsgLevel changeMsgLvl{RooFit::WARNING};
   RooWorkspace ws;
};

} // namespace

// The constraint terms of a pdf are the ones that depend on one of its parameters.
TEST_F(Constraints, ConstraintsOfPdf)
{
   RooArgSet all{pdf("cAlpha"), pdf("cBeta"), pdf("cMu"), pdf("cGamma")};
   RooArgSet observables{*ws.var("x"), *ws.cat("channel")};

   EXPECT_EQ(names(RooFit::constraintsOf(pdf("model"), all, &observables)), (Names{"cAlpha", "cBeta", "cMu"}));
   // Of a component of the model: only the ones of its own parameters.
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigA"), all, &observables)), (Names{"cAlpha", "cMu"}));
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigB"), all, &observables)), (Names{"cBeta", "cMu"}));
   // The observables can be left out if no term depends on them.
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigA"), all)), (Names{"cAlpha", "cMu"}));
   // Of a function.
   EXPECT_EQ(names(RooFit::constraintsOf(*ws.function("meanB"), all)), (Names{"cBeta", "cMu"}));

   // The term of a constant parameter is still a term of the model.
   ws.var("alpha")->setConstant(true);
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigA"), all, &observables)), (Names{"cAlpha", "cMu"}));

   // The order is the one of the constraint terms.
   RooArgList ordered{pdf("cMu"), pdf("cGamma"), pdf("cAlpha")};
   const RooArgSet selected = RooFit::constraintsOf(pdf("sigA"), ordered);
   ASSERT_EQ(selected.size(), 2u);
   EXPECT_STREQ(selected[0]->GetName(), "cMu");
   EXPECT_STREQ(selected[1]->GetName(), "cAlpha");

   // A term is no constraint term of itself, but another term of the same parameter is.
   ws.factory("Gaussian::cAlpha2(gAlpha2[0, -5, 5], alpha, 2.0)");
   all.add(pdf("cAlpha2"));
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("cAlpha"), all)), Names{"cAlpha2"});
}

// A pdf whose variable is an observable for one use and a parameter for another: the observables decide.
TEST_F(Constraints, ObservablesDecideWhatIsAParameter)
{
   ws.factory("Gaussian::cX(gX[0, -5, 5], x, 1.0)");
   RooArgSet all{pdf("cX"), pdf("cAlpha")};
   RooArgSet observables{*ws.var("x")};
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigA"), all, &observables)), Names{"cAlpha"});
   EXPECT_EQ(names(RooFit::constraintsOf(pdf("sigA"), all)), (Names{"cAlpha", "cX"}));
}

TEST_F(Constraints, ModelConfigsOfWorkspace)
{
   EXPECT_TRUE(RooFit::modelConfigs(ws).empty());

   addModelConfig("mc", "model");
   addModelConfig("mcA", "sigA");

   auto mcNames = [](std::vector<RooFit::ModelConfig *> const &mcs) {
      Names out;
      for (RooFit::ModelConfig const *mc : mcs) {
         out.insert(mc->GetName());
      }
      return out;
   };
   EXPECT_EQ(mcNames(RooFit::modelConfigs(ws)), (Names{"mc", "mcA"}));
   // The ones whose pdf is the given one, or has it as a component.
   EXPECT_EQ(mcNames(RooFit::modelConfigs(ws, &pdf("model"))), Names{"mc"});
   EXPECT_EQ(mcNames(RooFit::modelConfigs(ws, &pdf("sigA"))), (Names{"mc", "mcA"}));
   EXPECT_EQ(mcNames(RooFit::modelConfigs(ws, &pdf("sigB"))), Names{"mc"});
   EXPECT_EQ(mcNames(RooFit::modelConfigs(ws, ws.var("beta"))), Names{"mc"});
   EXPECT_TRUE(RooFit::modelConfigs(ws, &pdf("sigY")).empty());
}

// The constraint terms of pdfs that are part of no ModelConfig are registered in the workspace.
TEST_F(Constraints, RegisteredConstraintTerms)
{
   EXPECT_TRUE(RooFit::constraintTerms(ws).empty());
   EXPECT_TRUE(RooFit::constraintsOf(ws, pdf("model")).empty());

   RooFit::addConstraintTerms(ws, RooArgSet{pdf("cAlpha"), pdf("cGamma")});
   ASSERT_NE(ws.set(RooFit::constraintTermsSetName), nullptr);
   EXPECT_EQ(names(*ws.set(RooFit::constraintTermsSetName)), (Names{"cAlpha", "cGamma"}));

   // More terms extend the set, without a message about redefining it.
   {
      RooHelpers::HijackMessageStream hijack{RooFit::INFO, RooFit::InputArguments};
      RooFit::addConstraintTerms(ws, RooArgSet{pdf("cAlpha"), pdf("cBeta")});
      EXPECT_TRUE(hijack.str().empty()) << hijack.str();
   }
   EXPECT_EQ(names(RooFit::constraintTerms(ws)), (Names{"cAlpha", "cBeta", "cGamma"}));

   // A term that is not in the workspace is imported, and uses the parameters of the workspace.
   RooWorkspace other;
   other.factory("Gaussian::cOther(gOther[0, -5, 5], alpha[0, -5, 5], 1.0)");
   RooFit::addConstraintTerms(ws, RooArgSet{*other.pdf("cOther")});
   ASSERT_NE(ws.pdf("cOther"), nullptr);
   EXPECT_NE(ws.pdf("cOther"), other.pdf("cOther"));
   EXPECT_EQ(ws.pdf("cOther")->findServer("alpha"), ws.var("alpha"));

   // The constraint terms of a pdf are selected from the registered ones.
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("model"))), (Names{"cAlpha", "cBeta", "cOther"}));
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigB"))), Names{"cBeta"});
   EXPECT_TRUE(RooFit::constraintsOf(ws, pdf("sigY")).empty());
}

// A ModelConfig is what says which constraint terms its model has.
TEST_F(Constraints, ConstraintTermsOfModelConfigs)
{
   RooFit::addConstraintTerms(ws, RooArgSet{pdf("cAlpha"), pdf("cBeta")});
   const RooArgSet ofModel{pdf("cBeta"), pdf("cMu"), pdf("cGamma")};
   addModelConfig("mc", "model", &ofModel);

   // The model and its components: the external constraints of the ModelConfig, not the registered terms.
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("model"))), (Names{"cBeta", "cMu"}));
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigA"))), Names{"cMu"});
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigB"))), (Names{"cBeta", "cMu"}));
   // A pdf that is part of no ModelConfig: all terms that the workspace knows, the registered ones and the ones of
   // the ModelConfigs.
   ws.factory("Gaussian::sigAlpha(y, alpha, 1.0)");
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigAlpha"))), Names{"cAlpha"});
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigY"))), Names{"cMu"});

   // Several ModelConfigs with the pdf: all their terms.
   const RooArgSet ofChannel{pdf("cAlpha")};
   addModelConfig("mcA", "sigA", &ofChannel);
   EXPECT_EQ(names(RooFit::constraintsOf(ws, pdf("sigA"))), (Names{"cAlpha", "cMu"}));

   EXPECT_EQ(names(RooFit::constraintTerms(ws)), (Names{"cAlpha", "cBeta", "cMu", "cGamma"}));
}

TEST_F(Constraints, GenerateGlobalObservables)
{
   ws.var("alpha")->setVal(1.5);
   ws.var("beta")->setVal(-2.0);
   const RooArgSet constraints{pdf("cAlpha"), pdf("cBeta")};
   const RooArgSet globs{*ws.var("gAlpha"), *ws.var("gBeta"), *ws.var("gGamma")};

   const int n = 4000;
   std::unique_ptr<RooDataSet> data = RooFit::generateGlobalObservables(constraints, globs, n);
   ASSERT_NE(data, nullptr);
   EXPECT_EQ(data->numEntries(), n);
   // Only the global observables that the constraint terms depend on.
   EXPECT_EQ(names(*data->get()), (Names{"gAlpha", "gBeta"}));
   // Both terms have the width 1.
   EXPECT_NEAR(data->mean(*ws.var("gAlpha")), 1.5, 5. / std::sqrt(n));
   EXPECT_NEAR(data->mean(*ws.var("gBeta")), -2.0, 5. / std::sqrt(n));
   EXPECT_NEAR(data->sigma(*ws.var("gAlpha")), 1.0, 0.1);
   // The global observables themselves are not changed.
   EXPECT_DOUBLE_EQ(ws.var("gAlpha")->getVal(), 0.0);

   EXPECT_EQ(RooFit::generateGlobalObservables(constraints, globs)->numEntries(), 1);

   // Nothing to generate.
   EXPECT_EQ(RooFit::generateGlobalObservables(RooArgSet{}, globs), nullptr);
   EXPECT_EQ(RooFit::generateGlobalObservables(constraints, RooArgSet{}), nullptr);
   EXPECT_EQ(RooFit::generateGlobalObservables(constraints, RooArgSet{*ws.var("gGamma")}), nullptr);

   // The same for the sets of a ModelConfig.
   addModelConfig("mc", "model", &constraints);
   auto *mc = static_cast<RooFit::ModelConfig *>(ws.obj("mc"));
   EXPECT_EQ(mc->generateGlobalObservables(), nullptr);
   mc->SetGlobalObservables(globs);
   std::unique_ptr<RooDataSet> ofModel = mc->generateGlobalObservables(10);
   ASSERT_NE(ofModel, nullptr);
   EXPECT_EQ(ofModel->numEntries(), 10);
   EXPECT_EQ(names(*ofModel->get()), (Names{"gAlpha", "gBeta"}));
}

// A pdf without constraint terms stays what it is.
TEST_F(Constraints, SplitOffNoConstraints)
{
   RooArgSet observables{*ws.var("x"), *ws.var("y"), *ws.cat("channel")};
   for (const char *name : {"sigA", "model"}) {
      const RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf(name), observables);
      EXPECT_EQ(split.pdf, nullptr) << name;
      EXPECT_TRUE(split.constraints.empty()) << name;
   }
   // A product of pdfs of the observables.
   ws.factory("PROD::both(sigA, sigY)");
   const RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("both"), observables);
   EXPECT_EQ(split.pdf, nullptr);
   EXPECT_TRUE(split.constraints.empty());
}

// The constraint terms are the factors of a product that depend on no observable, but on a parameter of the other
// factors. The likelihood of what is left with the constraint terms is the one of the model.
TEST_F(Constraints, SplitOffConstraintsOfProduct)
{
   ws.factory("PROD::constrained(sigA, cAlpha, cMu)");
   RooRealVar &x = *ws.var("x");
   const RooArgSet observables{x};
   RooDataSet data{"data", "data", observables};
   for (double value : {-0.4, 0.3, 1.1}) {
      x.setVal(value);
      data.add(observables);
   }

   // The pdf with the constraint terms as factors can't be used for a likelihood, and the error names the way out.
   try {
      std::unique_ptr<RooAbsReal> nll{pdf("constrained").createNLL(data)};
      FAIL() << "expected an exception";
   } catch (std::runtime_error const &e) {
      EXPECT_NE(std::string(e.what()).find("RooFit::splitOffConstraints()"), std::string::npos) << e.what();
   }

   RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("constrained"), observables);
   EXPECT_EQ(names(split.constraints), (Names{"cAlpha", "cMu"}));
   // The constraint terms are the objects of the original pdf.
   EXPECT_TRUE(split.constraints.containsInstance(pdf("cAlpha")));
   // Only one factor is left of the product, so the pdf is that factor: a copy with the same parameters.
   ASSERT_NE(split.pdf, nullptr);
   EXPECT_STREQ(split.pdf->GetName(), "sigA");
   EXPECT_NE(split.pdf.get(), &pdf("sigA"));
   EXPECT_EQ(split.pdf->findServer("meanA"), ws.function("meanA"));
   // The original pdf is not changed.
   EXPECT_EQ(static_cast<RooProdPdf &>(pdf("constrained")).pdfList().size(), 3u);

   ws.var("alpha")->setVal(0.7);
   ws.var("mu")->setVal(-0.3);
   const RooArgSet globs{*ws.var("gAlpha"), *ws.var("gMu")};
   std::unique_ptr<RooAbsReal> nll{
      split.pdf->createNLL(data, RooFit::ExternalConstraints(split.constraints), RooFit::GlobalObservables(globs))};
   std::unique_ptr<RooAbsReal> nllOfData{pdf("sigA").createNLL(data)};
   double expected = nllOfData->getVal();
   for (const char *name : {"cAlpha", "cMu"}) {
      expected -= std::log(pdf(name).getVal(&globs));
   }
   EXPECT_NEAR(nll->getVal(), expected, 1e-12 * std::abs(expected));
}

// A factor is only a constraint term with respect to the observables.
TEST_F(Constraints, SplitOffDependsOnObservables)
{
   ws.factory("PROD::two(sigA, sigY, cAlpha)");

   // With both observables, two factors are left: a new product.
   {
      RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("two"), {*ws.var("x"), *ws.var("y")});
      EXPECT_EQ(names(split.constraints), Names{"cAlpha"});
      auto *prod = dynamic_cast<RooProdPdf *>(split.pdf.get());
      ASSERT_NE(prod, nullptr);
      EXPECT_STREQ(prod->GetName(), "two_unconstrained");
      EXPECT_EQ(names(prod->pdfList()), (Names{"sigA", "sigY"}));
      EXPECT_TRUE(prod->pdfList().containsInstance(pdf("sigA")));
      RooArgSet nset{*ws.var("x"), *ws.var("y")};
      EXPECT_DOUBLE_EQ(prod->getVal(&nset), pdf("sigA").getVal(&nset) * pdf("sigY").getVal(&nset));
   }
   // If only x is observed, the measurement of mu in y is a constraint term of mu.
   {
      RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("two"), {*ws.var("x")});
      EXPECT_EQ(names(split.constraints), (Names{"cAlpha", "sigY"}));
      ASSERT_NE(split.pdf, nullptr);
      EXPECT_STREQ(split.pdf->GetName(), "sigA");
   }
}

// A factor that shares no parameter with the others is no constraint term. The term of a constant parameter is one,
// although RooFit builds a likelihood for such a pdf.
TEST_F(Constraints, SplitOffSharedParameters)
{
   ws.factory("PROD::withOther(sigA, cAlpha, cGamma)");
   ws.var("alpha")->setConstant(true);
   const RooArgSet observables{*ws.var("x")};

   RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("withOther"), observables);
   EXPECT_EQ(names(split.constraints), Names{"cAlpha"});
   auto *prod = dynamic_cast<RooProdPdf *>(split.pdf.get());
   ASSERT_NE(prod, nullptr);
   EXPECT_EQ(names(prod->pdfList()), (Names{"sigA", "cGamma"}));

   RooDataSet data{"data", "data", observables};
   data.add(observables);
   EXPECT_NO_THROW(std::unique_ptr<RooAbsReal>{pdf("withOther").createNLL(data)});
   EXPECT_NO_THROW(std::unique_ptr<RooAbsReal>{split.pdf->createNLL(data)});
}

// The channels of a RooSimultaneous, with a constraint term that they share and a product that groups terms.
TEST_F(Constraints, SplitOffConstraintsOfSimultaneous)
{
   ws.factory("PROD::chA(sigA, cAlpha, cMu)");
   ws.factory("PROD::chB(sigB, PROD::grouped(cBeta, cMu))");
   ws.factory("SIMUL::legacy(channel, A=chA, B=chB)");
   const RooArgSet observables{*ws.var("x"), *ws.cat("channel")};

   RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("legacy"), observables);
   // Every term once, and the terms of the group instead of the group.
   EXPECT_EQ(names(split.constraints), (Names{"cAlpha", "cBeta", "cMu"}));
   auto *sim = dynamic_cast<RooSimultaneous *>(split.pdf.get());
   ASSERT_NE(sim, nullptr);
   EXPECT_STREQ(sim->GetName(), "legacy_unconstrained");
   EXPECT_EQ(&sim->indexCat(), ws.cat("channel"));
   // The channels are the pdfs of the original model themselves.
   EXPECT_EQ(sim->getPdf("A"), &pdf("sigA"));
   EXPECT_EQ(sim->getPdf("B"), &pdf("sigB"));
   // The original model is not changed.
   EXPECT_EQ(static_cast<RooSimultaneous &>(pdf("legacy")).getPdf("A"), &pdf("chA"));
   EXPECT_EQ(static_cast<RooProdPdf &>(pdf("chA")).pdfList().size(), 3u);

   // The same likelihood as the model that is built without constraint terms in the pdf.
   std::unique_ptr<RooDataSet> data = makeData();
   EXPECT_THROW(std::unique_ptr<RooAbsReal>{pdf("legacy").createNLL(*data)}, std::runtime_error);
   const RooArgSet globs{*ws.var("gAlpha"), *ws.var("gBeta"), *ws.var("gMu")};
   const RooArgSet constraints{pdf("cAlpha"), pdf("cBeta"), pdf("cMu")};
   std::unique_ptr<RooAbsReal> nll{
      split.pdf->createNLL(*data, RooFit::ExternalConstraints(split.constraints), RooFit::GlobalObservables(globs))};
   std::unique_ptr<RooAbsReal> expected{
      pdf("model").createNLL(*data, RooFit::ExternalConstraints(constraints), RooFit::GlobalObservables(globs))};
   for (auto const &[alpha, beta, mu] :
        std::vector<std::tuple<double, double, double>>{{0., 0., 0.}, {0.7, -1.2, 0.4}, {-2.0, 0.3, 1.1}}) {
      ws.var("alpha")->setVal(alpha);
      ws.var("beta")->setVal(beta);
      ws.var("mu")->setVal(mu);
      EXPECT_NEAR(nll->getVal(), expected->getVal(), 1e-12 * std::abs(expected->getVal())) << alpha;
   }
}

// A channel that is used for several states is replaced by one pdf.
TEST_F(Constraints, SplitOffConstraintsOfSharedChannel)
{
   ws.factory("PROD::both(sigA, sigY, cAlpha)");
   ws.factory("SIMUL::shared(region[R1, R2], R1=both, R2=both)");
   RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("shared"), {*ws.var("x"), *ws.var("y"), *ws.cat("region")});
   EXPECT_EQ(names(split.constraints), Names{"cAlpha"});
   auto *sim = dynamic_cast<RooSimultaneous *>(split.pdf.get());
   ASSERT_NE(sim, nullptr);
   ASSERT_NE(sim->getPdf("R1"), nullptr);
   EXPECT_STREQ(sim->getPdf("R1")->GetName(), "both_unconstrained");
   EXPECT_EQ(sim->getPdf("R1"), sim->getPdf("R2"));
}

// Products in products, and the pdf of a RooExtendPdf.
TEST_F(Constraints, SplitOffConstraintsOfNestedPdfs)
{
   ws.factory("PROD::inner(sigA, sigY, cAlpha)");
   ws.factory("PROD::outer(inner, cMu)");
   const RooArgSet observables{*ws.var("x"), *ws.var("y")};
   {
      RooFit::SplitPdf split = RooFit::splitOffConstraints(pdf("outer"), observables);
      EXPECT_EQ(names(split.constraints), (Names{"cAlpha", "cMu"}));
      auto *prod = dynamic_cast<RooProdPdf *>(split.pdf.get());
      ASSERT_NE(prod, nullptr);
      EXPECT_STREQ(prod->GetName(), "inner_unconstrained");
      EXPECT_EQ(names(prod->pdfList()), (Names{"sigA", "sigY"}));
   }

   RooRealVar n{"n", "n", 100., 0., 1000.};
   RooExtendPdf extended{"extended", "extended", pdf("outer"), n};
   {
      RooFit::SplitPdf split = RooFit::splitOffConstraints(extended, observables);
      EXPECT_EQ(names(split.constraints), (Names{"cAlpha", "cMu"}));
      auto *ext = dynamic_cast<RooExtendPdf *>(split.pdf.get());
      ASSERT_NE(ext, nullptr);
      EXPECT_STREQ(ext->GetName(), "extended_unconstrained");
      EXPECT_STREQ(ext->pdf().GetName(), "inner_unconstrained");
      EXPECT_DOUBLE_EQ(ext->expectedEvents(&observables), 100.);
      // The new pdf owns the nodes that were created for it: it can be used on its own.
      RooArgSet nset{observables};
      EXPECT_DOUBLE_EQ(ext->getVal(&nset), pdf("sigA").getVal(&nset) * pdf("sigY").getVal(&nset));
   }
}

// A product can make its factors conditional pdfs. This is kept, also if only one factor is left.
TEST_F(Constraints, SplitOffKeepsConditionalFactors)
{
   RooProdPdf conditional{"conditional", "conditional", RooArgSet{pdf("cAlpha")},
                          RooFit::Conditional(RooArgSet{pdf("sigA")}, RooArgSet{*ws.var("x")})};
   RooFit::SplitPdf split = RooFit::splitOffConstraints(conditional, {*ws.var("x")});
   EXPECT_EQ(names(split.constraints), Names{"cAlpha"});
   auto *prod = dynamic_cast<RooProdPdf *>(split.pdf.get());
   ASSERT_NE(prod, nullptr);
   EXPECT_STREQ(prod->GetName(), "conditional_unconstrained");
   EXPECT_EQ(names(prod->pdfList()), Names{"sigA"});
   RooArgSet *nset = prod->findPdfNSet(pdf("sigA"));
   ASSERT_NE(nset, nullptr);
   EXPECT_EQ(names(*nset), Names{"x"});
}
