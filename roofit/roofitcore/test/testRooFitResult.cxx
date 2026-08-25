// Tests for the RooFitResult
// Authors: Jonas Rembser, CERN 2026

#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHelpers.h>
#include <RooPolyVar.h>
#include <RooRandom.h>
#include <RooRealVar.h>

#include <TMatrixDSym.h>

#include <gtest/gtest.h>

#include <memory>

namespace {

/// Create a fit result that represents the measurement of two parameters
/// "p0" and "p1" with central values (10, 20) and covariance {{4, 1}, {1, 9}}.
std::unique_ptr<RooFitResult> makeTestFitResult(RooArgList const &params)
{
   auto fr = std::make_unique<RooFitResult>("fr", "fr");
   fr->setConstParList({});
   fr->setInitParList(params);
   fr->setFinalParList(params);
   TMatrixDSym cov(2);
   cov(0, 0) = 4.0;
   cov(0, 1) = 1.0;
   cov(1, 0) = 1.0;
   cov(1, 1) = 9.0;
   TMatrixDSym corr(2);
   corr(0, 0) = 1.0;
   corr(0, 1) = cov(0, 1) / std::sqrt(cov(0, 0) * cov(1, 1));
   corr(1, 0) = corr(0, 1);
   corr(1, 1) = 1.0;
   fr->fillCorrMatrix({corr(0, 1), corr(0, 1)}, corr, cov);
   fr->setStatus(0);
   fr->setCovQual(3);
   return fr;
}

/// The chi-square corresponding to the fit result from makeTestFitResult()
/// with the predictions {theta, 2 * theta}, computed with the explicit
/// inverse covariance matrix.
double referenceChi2(double theta)
{
   const double d0 = 10. - theta;
   const double d1 = 20. - 2. * theta;
   // inverse of {{4, 1}, {1, 9}} is {{9, -1}, {-1, 4}} / 35
   return (9. * d0 * d0 - 2. * d0 * d1 + 4. * d1 * d1) / 35.;
}

} // namespace

/// Check that the chi-square model created from a fit result reproduces the
/// analytically-known weighted least squares solution for linear predictions.
TEST(RooFitResult, CreateChi2Pdf)
{
   RooRealVar p0{"p0", "p0", 10., 0., 100.};
   RooRealVar p1{"p1", "p1", 20., 0., 100.};
   RooArgList params{p0, p1};

   std::unique_ptr<RooFitResult> fr = makeTestFitResult(params);

   // Predictions linear in theta: mu0 = theta, mu1 = 2 * theta
   RooRealVar theta{"theta", "theta", 1., -100., 100.};
   RooPolyVar pred1{"pred1", "pred1", theta, RooArgList{0.0, 2.0}};

   std::unique_ptr<RooAbsPdf> chi2Pdf = fr->createChi2Pdf(params, {theta, pred1});
   ASSERT_NE(chi2Pdf, nullptr);

   std::unique_ptr<RooDataSet> data = fr->createChi2DataSet(params);
   ASSERT_NE(data, nullptr);
   ASSERT_EQ(data->numEntries(), 1);
   RooArgSet const &row = *data->get(0);
   EXPECT_DOUBLE_EQ(static_cast<RooAbsReal const &>(row["p0_obs"]).getVal(), 10.);
   EXPECT_DOUBLE_EQ(static_cast<RooAbsReal const &>(row["p1_obs"]).getVal(), 20.);

   // Differences of twice the negative log-likelihood must match the
   // analytically computed chi-square differences.
   std::unique_ptr<RooAbsReal> nll{chi2Pdf->createNLL(*data)};
   auto nllVal = [&](double th) {
      theta.setVal(th);
      return nll->getVal();
   };
   const double nllRef = nllVal(10.);
   for (double th : {5., 7.5, 12., 20.}) {
      EXPECT_NEAR(2. * (nllVal(th) - nllRef), referenceChi2(th) - referenceChi2(10.), 1e-10);
   }

   // The analytic weighted least squares solution for this linear model is
   // exactly theta = 10 with variance (A^T V^-1 A)^-1 = 5/3. The tolerance
   // on the fitted value reflects Minuit's default EDM stopping condition,
   // which is not very strict relative to the large parameter uncertainty of
   // about 1.29.
   theta.setVal(1.);
   theta.setError(1.);
   std::unique_ptr<RooFitResult> frChi2{chi2Pdf->fitTo(*data, RooFit::Save(), RooFit::PrintLevel(-1))};
   ASSERT_EQ(frChi2->status(), 0);
   EXPECT_NEAR(theta.getVal(), 10., 0.05);
   EXPECT_NEAR(theta.getError(), std::sqrt(5. / 3.), 1e-3 * theta.getError());

   // The result must not depend on the order in which the parameters are
   // given, as long as the predictions are matched by position.
   std::unique_ptr<RooAbsPdf> chi2PdfRev = fr->createChi2Pdf({p1, p0}, {pred1, theta});
   ASSERT_NE(chi2PdfRev, nullptr);
   std::unique_ptr<RooDataSet> dataRev = fr->createChi2DataSet({p1, p0});
   std::unique_ptr<RooAbsReal> nllRev{chi2PdfRev->createNLL(*dataRev)};
   for (double th : {5., 10., 20.}) {
      theta.setVal(th);
      EXPECT_NEAR(nllRev->getVal(), nll->getVal(), 1e-10);
   }
}

/// Check that toys generated from the chi-square model follow the covariance
/// matrix stored in the fit result.
TEST(RooFitResult, CreateChi2PdfGenerate)
{
   RooRealVar p0{"p0", "p0", 10., 0., 100.};
   RooRealVar p1{"p1", "p1", 20., 0., 100.};
   RooArgList params{p0, p1};

   std::unique_ptr<RooFitResult> fr = makeTestFitResult(params);

   RooRealVar theta{"theta", "theta", 10., -100., 100.};
   RooPolyVar pred1{"pred1", "pred1", theta, RooArgList{0.0, 2.0}};

   std::unique_ptr<RooAbsPdf> chi2Pdf = fr->createChi2Pdf(params, {theta, pred1});
   ASSERT_NE(chi2Pdf, nullptr);
   std::unique_ptr<RooDataSet> data = fr->createChi2DataSet(params);
   std::unique_ptr<RooArgSet> obs{chi2Pdf->getObservables(*data)};
   ASSERT_EQ(obs->size(), 2);

   RooRandom::randomGenerator()->SetSeed(1337);
   const std::size_t nToys = 10000;
   std::unique_ptr<RooDataSet> toys{chi2Pdf->generate(*obs, nToys)};
   ASSERT_EQ(toys->numEntries(), int(nToys));

   // At theta = 10, the generation means are mu = (10, 20).
   double mean[2] = {0., 0.};
   double cov[3] = {0., 0., 0.}; // V(0,0), V(0,1), V(1,1)
   for (int i = 0; i < toys->numEntries(); ++i) {
      RooArgSet const &row = *toys->get(i);
      mean[0] += static_cast<RooAbsReal const &>(row["p0_obs"]).getVal();
      mean[1] += static_cast<RooAbsReal const &>(row["p1_obs"]).getVal();
   }
   mean[0] /= nToys;
   mean[1] /= nToys;
   for (int i = 0; i < toys->numEntries(); ++i) {
      RooArgSet const &row = *toys->get(i);
      const double d0 = static_cast<RooAbsReal const &>(row["p0_obs"]).getVal() - mean[0];
      const double d1 = static_cast<RooAbsReal const &>(row["p1_obs"]).getVal() - mean[1];
      cov[0] += d0 * d0;
      cov[1] += d0 * d1;
      cov[2] += d1 * d1;
   }
   for (double &element : cov) {
      element /= nToys - 1;
   }

   EXPECT_NEAR(mean[0], 10., 5. * 2. / std::sqrt(nToys));
   EXPECT_NEAR(mean[1], 20., 5. * 3. / std::sqrt(nToys));
   EXPECT_NEAR(cov[0], 4., 0.2);
   EXPECT_NEAR(cov[1], 1., 0.2);
   EXPECT_NEAR(cov[2], 9., 0.5);
}

/// Check that invalid inputs are rejected.
TEST(RooFitResult, CreateChi2PdfInputValidation)
{
   RooHelpers::LocalChangeMsgLevel silence{RooFit::FATAL};

   RooRealVar p0{"p0", "p0", 10., 0., 100.};
   RooRealVar p1{"p1", "p1", 20., 0., 100.};
   RooArgList params{p0, p1};

   std::unique_ptr<RooFitResult> fr = makeTestFitResult(params);

   RooRealVar theta{"theta", "theta", 1., -100., 100.};
   RooRealVar other{"other", "other", 1., 0., 2.};

   // Mismatch between number of parameters and predictions
   EXPECT_EQ(fr->createChi2Pdf(params, {theta}), nullptr);
   // Empty inputs
   EXPECT_EQ(fr->createChi2Pdf({}, {}), nullptr);
   EXPECT_EQ(fr->createChi2DataSet({}), nullptr);
   // Parameter that was not floating in the fit
   EXPECT_EQ(fr->createChi2Pdf({p0, other}, {theta, theta}), nullptr);
   EXPECT_EQ(fr->createChi2DataSet({p0, other}), nullptr);
   // Duplicate parameter
   EXPECT_EQ(fr->createChi2Pdf({p0, p0}, {theta, theta}), nullptr);
}
