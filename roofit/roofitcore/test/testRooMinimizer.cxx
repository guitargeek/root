// Tests for the RooProdPdf
// Author: Jonas Rembser, CERN, October 2024

#include <RooAddition.h>
#include <RooAddPdf.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooHelpers.h>
#include <RooMinimizer.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>

#include <TMatrixDSym.h>

#include <cmath>

#include "gtest_wrapper.h"

class EvalBackendParametrizedTest : public testing::TestWithParam<std::tuple<RooFit::EvalBackend>> {
public:
   EvalBackendParametrizedTest() : _evalBackend{RooFit::EvalBackend(RooFit::EvalBackend::Value::Legacy)} {}

private:
   void SetUp() override
   {
      RooRandom::randomGenerator()->SetSeed(1337ul);
      _evalBackend = std::get<0>(GetParam());
      _changeMsgLvl = std::make_unique<RooHelpers::LocalChangeMsgLevel>(RooFit::WARNING);
   }

   void TearDown() override { _changeMsgLvl.reset(); }

protected:
   RooFit::EvalBackend _evalBackend;

private:
   std::unique_ptr<RooHelpers::LocalChangeMsgLevel> _changeMsgLvl;
};

// Unit test corresponding to the rf601 tutorial, but parametrized for the
// different evaluation backends.
TEST_P(EvalBackendParametrizedTest, RF601)
{
   RooRealVar x("x", "x", -20, 20);

   // Model (intentional strong correlations)
   RooRealVar mean("mean", "mean of g1 and g2", 0);
   RooRealVar sigma_g1("sigma_g1", "width of g1", 3, 1.0, 5.0);
   RooGaussian g1("g1", "g1", x, mean, sigma_g1);

   RooRealVar sigma_g2("sigma_g2", "width of g2", 4, 3.0, 6.0);
   RooGaussian g2("g2", "g2", x, mean, sigma_g2);

   RooRealVar frac("frac", "frac", 0.5, 0.0, 1.0);
   RooAddPdf model("model", "model", RooArgList(g1, g2), frac);

   std::unique_ptr<RooDataSet> data{model.generate(x, 1000)};

   std::unique_ptr<RooAbsReal> nll{model.createNLL(*data, RooFit::EvalBackend(_evalBackend))};

   // Reference fit results. We are building them manually in this code in
   // order to avoid binary reference files.
   RooFitResult rRef;
   RooFitResult r2Ref;

   rRef.setInitParList({frac, sigma_g1, sigma_g2});
   r2Ref.setInitParList({frac, sigma_g1});

   {
      RooArgList params{frac, sigma_g1, sigma_g2};
      rRef.setMinNLL(2690.1655479975625);
      std::vector<double> valsRef{0.48732473752392391, 2.7579472102402196, 4.2132158568489011};
      std::vector<double> errorsRef{0.29415767363636769, 0.48854402667226893, 0.47322653666757764};
      std::vector<double> globalCC{0.98208607303266726, 0.94654805382506446, 0.9542122909792754};
      std::vector<double> corrsV{1.0, 0.927913, 0.938339, 0.927913, 1.0, 0.806083, 0.938339, 0.806083, 1.0};
      const std::size_t nParams = params.size();

      RooArgList paramsPostFit;
      params.snapshot(paramsPostFit, false);

      for (std::size_t i = 0; i < nParams; ++i) {
         auto &var = static_cast<RooRealVar &>(paramsPostFit[i]);
         var.setVal(valsRef[i]);
         var.setError(errorsRef[i]);
      }

      rRef.setConstParList(mean);
      rRef.setFinalParList(paramsPostFit);
      TMatrixDSym corrs(nParams);
      TMatrixDSym covs(nParams);
      corrs.SetMatrixArray(corrsV.data());
      rRef.fillCorrMatrix(globalCC, corrs, covs);
   }

   RooMinimizer m(*nll);

   m.setPrintLevel(-1);

   m.migrad();
   m.hesse();

   m.minos(sigma_g2); // Run MINOS on sigma_g2 parameter only

   std::unique_ptr<RooFitResult> r{m.save()};

   // You can manually change the value of a (constant) parameter
   mean = 0.3;

   {
      RooArgList params{frac, sigma_g1};
      r2Ref.setMinNLL(2698.6818728208696);
      std::vector<double> constValsRef{0.3, 4.1869317966930382};
      std::vector<double> constErrorsRef{0.0, 0.0};
      std::vector<double> valsRef{0.45212322018672196, 2.746677666046311};
      std::vector<double> errorsRef{0.11175880034776256, 0.31839226266447529};
      std::vector<double> globalCC{0.84401375662480127, 0.84401375662480116, 0.};
      // clang-format off
      std::vector<double> corrsV{1., 0.84401326379874853, 0,
                                 0.84401326379874853, 1., 0.,
                                 0., 0., 1.};
      // clang-format on

      RooArgList constParams;
      RooArgList{mean, sigma_g2}.snapshot(constParams, false);

      for (std::size_t i = 0; i < constParams.size(); ++i) {
         auto &var = static_cast<RooRealVar &>(constParams[i]);
         var.setVal(constValsRef[i]);
         var.setError(constErrorsRef[i]);
      }

      RooArgList paramsPostFit;
      params.snapshot(paramsPostFit, false);

      for (std::size_t i = 0; i < paramsPostFit.size(); ++i) {
         auto &var = static_cast<RooRealVar &>(paramsPostFit[i]);
         var.setVal(valsRef[i]);
         var.setError(errorsRef[i]);
      }

      r2Ref.setConstParList(constParams);
      r2Ref.setFinalParList(paramsPostFit);
      TMatrixDSym corrs(3);
      TMatrixDSym covs(3);
      corrs.SetMatrixArray(corrsV.data());
      r2Ref.fillCorrMatrix(globalCC, corrs, covs);
   }

   // Rerun MIGRAD,HESSE
   m.migrad();
   m.hesse();

   // Now fix sigma_g2
   sigma_g2.setConstant(true);

   // Rerun MIGRAD,HESSE
   m.migrad();
   m.hesse();

   std::unique_ptr<RooFitResult> r2{m.save()};

   // The tolerance parameter is necessary because not all backends give
   // exactly the same results: when using AD, the final result is slightly
   // different.
   const double tol = 1e-4;

   EXPECT_TRUE(r->isIdentical(rRef, tol, tol));
   EXPECT_TRUE(r2->isIdentical(r2Ref, tol, tol));
}

INSTANTIATE_TEST_SUITE_P(RooMinimizer, EvalBackendParametrizedTest, testing::Values(ROOFIT_EVAL_BACKENDS_WITH_CODEGEN),
                         [](testing::TestParamInfo<EvalBackendParametrizedTest::ParamType> const &paramInfo) {
                            std::stringstream ss;
                            ss << "EvalBackend" << std::get<0>(paramInfo.param).name();
                            return ss.str();
                         });

// Check the vanishing-second-derivative optimization in MnHesse, which is
// driven by the parameter independence information that RooMinimizerFcn
// derives from the computation graph. The covariance matrix from hesse() must
// agree with the analytical one, both when all parameters float and when a
// parameter is fixed after the RooMinimizer was constructed. The latter case
// is a regression test for the translation between Minuit-internal parameter
// indices (which exclude fixed parameters) and the external indices that
// RooMinimizerFcn::secondDerivativeAlwaysVanishes() is defined in.
TEST(RooMinimizer, SecondDerivativeAlwaysVanishesHesse)
{
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};

   RooRealVar a("a", "a", 1, -10, 10);
   RooRealVar b("b", "b", 1, -10, 10);
   RooRealVar c("c", "c", 1, -10, 10);

   // Two additive terms: t1 couples (a,c) and t2 couples (b,c), so only the
   // mixed second derivative w.r.t. (a,b) vanishes identically and can be
   // skipped in the numerical Hessian.
   RooFormulaVar t1("t1", "t1", "(a-1)^2 + (a-c)^2 + (c-1)^2", {a, c});
   RooFormulaVar t2("t2", "t2", "(b-c)^2 + (b-1)^2", {b, c});
   RooAddition f("f", "f", {t1, t2});

   auto resetParams = [&]() {
      for (RooRealVar *var : {&a, &b, &c}) {
         var->setVal(0.5);
         var->setError(0.0);
         var->setConstant(false);
      }
   };

   // The function is quadratic, so the finite-difference Hessian is exact up
   // to numerical noise.
   const double tol = 1e-3;

   // All parameters floating: with errorDef = 1, the covariance is 2 * H^-1
   // with H = [[4, 0, -2], [0, 4, -2], [-2, -2, 6]] in the order (a, b, c).
   {
      resetParams();
      RooMinimizer m(f);
      m.setPrintLevel(-1);
      m.migrad();
      m.hesse();
      EXPECT_NEAR(a.getError(), std::sqrt(0.625), tol);
      EXPECT_NEAR(b.getError(), std::sqrt(0.625), tol);
      EXPECT_NEAR(c.getError(), std::sqrt(0.5), tol);
   }

   // With "a" fixed, the reduced Hessian in (b, c) is [[4, -2], [-2, 6]],
   // so the covariance is 2 * H^-1 = [[0.6, 0.2], [0.2, 0.4]]. This must not
   // depend on whether "a" was fixed before or after constructing the
   // minimizer: in the latter case, Minuit's internal parameter indices no
   // longer coincide with the external ones that the vanishing-second-
   // derivative mask is defined in, so wrong bookkeeping would zero out the
   // genuinely non-vanishing (b, c) element here.
   for (bool fixAfterConstruction : {false, true}) {
      resetParams();
      a.setVal(1.0);
      if (!fixAfterConstruction) {
         a.setConstant(true);
      }
      RooMinimizer m(f);
      m.setPrintLevel(-1);
      if (fixAfterConstruction) {
         a.setConstant(true);
      }
      m.migrad();
      m.hesse();
      EXPECT_NEAR(b.getError(), std::sqrt(0.6), tol) << "fixAfterConstruction = " << fixAfterConstruction;
      EXPECT_NEAR(c.getError(), std::sqrt(0.4), tol) << "fixAfterConstruction = " << fixAfterConstruction;
   }
}

// The same optimization for a simultaneous fit, whose likelihood
// RooSimultaneous::compileForNormSet() compiles into a mixture over the
// channels. The per-channel additive structure is not a RooAddition there but
// the gated terms of that mixture, and missing it silently gives up all
// cross-channel independence, which is exactly the many-channel case the
// optimization is meant for.
//
// Both directions matter and are covered here: independence that is there must
// be exploited, and coupling that is there must survive. The references are
// analytical, so they are independent of RooFit's own Hessian.
TEST(RooMinimizer, SecondDerivativeAlwaysVanishesHesseSimultaneous)
{
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING};
   RooRandom::randomGenerator()->SetSeed(1337ul);

   const int nEvents = 5000;

   // Part 1: channels that share no parameter at all. For an (extended)
   // Gaussian fit, the observed information at the maximum-likelihood point is
   // exact: var(mean) = s^2/N, var(sigma) = s^2/(2N) and var(yield) = N, with
   // vanishing mixed terms. The range is wide enough that the truncation of
   // the Gaussian normalization is far below the tolerance.
   {
      RooRealVar x("x", "x", -30, 30);
      RooCategory cat("cat", "cat");
      RooSimultaneous sim("sim", "sim", cat);

      std::vector<std::unique_ptr<RooAbsArg>> keep;
      std::vector<RooRealVar *> means;
      std::vector<RooRealVar *> sigmas;
      std::vector<RooRealVar *> yields;

      RooDataSet data{"data", "data", {x, cat}};

      for (int i = 0; i < 3; ++i) {
         const std::string ch = "ch" + std::to_string(i);
         cat.defineType(ch.c_str(), i);

         auto mean = std::make_unique<RooRealVar>(("m_" + ch).c_str(), "", 0.5 * i, -10, 10);
         auto sigma = std::make_unique<RooRealVar>(("s_" + ch).c_str(), "", 1.0 + 0.2 * i, 0.1, 10);
         auto gauss = std::make_unique<RooGaussian>(("g_" + ch).c_str(), "", x, *mean, *sigma);
         auto yield = std::make_unique<RooRealVar>(("n_" + ch).c_str(), "", nEvents, 0., 10. * nEvents);
         auto ext = std::make_unique<RooExtendPdf>(("e_" + ch).c_str(), "", *gauss, *yield);

         cat.setIndex(i);
         std::unique_ptr<RooDataSet> chanData{gauss->generate(x, nEvents)};
         for (int j = 0; j < chanData->numEntries(); ++j) {
            x.setVal(chanData->get(j)->getRealValue("x"));
            data.add({x, cat});
         }

         sim.addPdf(*ext, ch.c_str());
         means.push_back(mean.get());
         sigmas.push_back(sigma.get());
         yields.push_back(yield.get());
         keep.emplace_back(std::move(mean));
         keep.emplace_back(std::move(sigma));
         keep.emplace_back(std::move(gauss));
         keep.emplace_back(std::move(yield));
         keep.emplace_back(std::move(ext));
      }

      std::unique_ptr<RooAbsReal> nll{sim.createNLL(data)};
      RooMinimizer m{*nll};
      m.setPrintLevel(-1);
      m.migrad();
      m.zeroEvalCount();
      m.hesse();
      const int nEvals = m.evalCounter();

      for (std::size_t i = 0; i < means.size(); ++i) {
         const double s = sigmas[i]->getVal();
         const double n = yields[i]->getVal();
         EXPECT_NEAR(means[i]->getError(), s / std::sqrt(n), 1e-3 * s / std::sqrt(n)) << "channel " << i;
         EXPECT_NEAR(sigmas[i]->getError(), s / std::sqrt(2. * n), 1e-3 * s / std::sqrt(2. * n)) << "channel " << i;
         EXPECT_NEAR(yields[i]->getError(), std::sqrt(n), 1e-3 * std::sqrt(n)) << "channel " << i;
      }

      // The values above are also what the unoptimized Hessian gives, so check
      // that the independence was actually exploited. MnHesse spends one
      // function evaluation per parameter pair it does not skip; of the 36
      // pairs of the nine parameters here, the 27 that cross a channel
      // boundary are skipped, which takes the count from 73 down to 46. The
      // bound leaves room for MnHesse's own bookkeeping to change.
      EXPECT_LT(nEvals, 60);
   }

   // Part 2: a parameter shared by all channels must stay coupled to each of
   // them. With the widths held constant, the channel means m_i + shift make
   // the likelihood an exact quadratic form in (shift, m_1, m_2), so the
   // covariance is the inverse of a Hessian that can be written down: with
   // a_i = N_i / s_i^2, the mean of the anchor channel held constant, and the
   // parameters ordered (shift, m_1, m_2),
   //
   //     H = [[a_0 + a_1 + a_2, a_1, a_2], [a_1, a_1, 0], [a_2, 0, a_2]].
   //
   // The off-diagonal (shift, m_i) entries are large, so wrongly advertising
   // them as vanishing would show up immediately.
   {
      RooRealVar x("x", "x", -30, 30);
      RooCategory cat("cat", "cat");
      RooSimultaneous sim("sim", "sim", cat);
      RooRealVar shift("shift", "", 0.0, -5, 5);

      std::vector<std::unique_ptr<RooAbsArg>> keep;
      std::vector<RooRealVar *> offsets;
      std::vector<double> a;

      RooDataSet data{"data", "data", {x, cat}};

      for (int i = 0; i < 3; ++i) {
         const std::string ch = "ch" + std::to_string(i);
         cat.defineType(ch.c_str(), i);

         auto offset = std::make_unique<RooRealVar>(("m_" + ch).c_str(), "", 0.4 * i, -10, 10);
         // The first channel anchors the absolute scale, so that the common
         // shift is identifiable instead of exactly degenerate with the means.
         offset->setConstant(i == 0);
         auto mean = std::make_unique<RooFormulaVar>(("mean_" + ch).c_str(), "", "@0+@1", RooArgList{*offset, shift});
         auto sigma = std::make_unique<RooRealVar>(("s_" + ch).c_str(), "", 1.0 + 0.2 * i, 0.1, 10);
         sigma->setConstant(true);
         auto gauss = std::make_unique<RooGaussian>(("g_" + ch).c_str(), "", x, *mean, *sigma);

         cat.setIndex(i);
         std::unique_ptr<RooDataSet> chanData{gauss->generate(x, nEvents)};
         for (int j = 0; j < chanData->numEntries(); ++j) {
            x.setVal(chanData->get(j)->getRealValue("x"));
            data.add({x, cat});
         }

         sim.addPdf(*gauss, ch.c_str());
         a.push_back(nEvents / (sigma->getVal() * sigma->getVal()));
         if (i != 0) {
            offsets.push_back(offset.get());
         }
         keep.emplace_back(std::move(offset));
         keep.emplace_back(std::move(mean));
         keep.emplace_back(std::move(sigma));
         keep.emplace_back(std::move(gauss));
      }

      std::unique_ptr<RooAbsReal> nll{sim.createNLL(data)};
      RooMinimizer m{*nll};
      m.setPrintLevel(-1);
      m.migrad();
      m.hesse();
      std::unique_ptr<RooFitResult> res{m.save()};

      TMatrixDSym hessian(3);
      hessian(0, 0) = a[0] + a[1] + a[2];
      hessian(0, 1) = hessian(1, 0) = a[1];
      hessian(0, 2) = hessian(2, 0) = a[2];
      hessian(1, 1) = a[1];
      hessian(2, 2) = a[2];
      hessian(1, 2) = hessian(2, 1) = 0.0;
      TMatrixDSym reference{hessian};
      reference.Invert();

      EXPECT_NEAR(shift.getError(), std::sqrt(reference(0, 0)), 1e-3 * std::sqrt(reference(0, 0)));
      for (std::size_t i = 0; i < offsets.size(); ++i) {
         const double ref = std::sqrt(reference(i + 1, i + 1));
         EXPECT_NEAR(offsets[i]->getError(), ref, 1e-3 * ref) << "channel " << i + 1;
      }

      // The correlation that couples the channels must not be lost, while the
      // two channel means stay uncorrelated with each other.
      const RooArgList &pars = res->floatParsFinal();
      auto corr = [&](const char *n1, const char *n2) {
         return res->correlation(pars.find(n1)->GetName(), pars.find(n2)->GetName());
      };
      const double refCorr01 = reference(0, 1) / std::sqrt(reference(0, 0) * reference(1, 1));
      EXPECT_NEAR(corr("shift", "m_ch1"), refCorr01, 1e-3);
      EXPECT_NEAR(corr("m_ch1", "m_ch2"), reference(1, 2) / std::sqrt(reference(1, 1) * reference(2, 2)), 1e-3);
   }
}
