#include <RooCategory.h>
#include <RooConstVar.h>
#include <RooGaussian.h>
#include <RooMultiPdf.h>
#include <RooMultiReal.h>
#include <RooRealVar.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
#include <RooDataSet.h>
#include <RooRandom.h>
#include <RooMinimizer.h>
#include <gtest/gtest.h>

#include <memory>

// Helper function to count parameters including x
int countFloatingParametersIncludingObservable(const RooAbsPdf &pdf)
{
   std::unique_ptr<RooArgSet> params{pdf.getParameters(RooArgSet())};
   int count = 0;
   for (auto *obj : *params) {
      auto *var = dynamic_cast<RooRealVar *>(obj);
      if (var && !var->isConstant()) {
         count++;
      }
   }

   return count;
}

TEST(RooMultiPdf, SelectsCorrectPdf)
{

   RooRealVar x("x", "x", -10, 10);
   x.setVal(2.0);

   RooRealVar m1("mean1", "mean1", 0.);
   RooRealVar s1("sigma1", "sigma1", 1., 0.001, 10.);
   RooRealVar m2("mean2", "mean2", 2.);
   RooRealVar s2("sigma2", "sigma2", 1., 0.001, 10.);
   RooGaussian gaus1("gaus1", "gaus1", x, m1, s1);
   RooGaussian gaus2("gaus2", "gaus2", x, m2, s2);

   RooCategory indx("my_special_index", "my_index");

   RooArgList list{gaus1, gaus2};

   RooMultiPdf pdf("mult", "multi_pdf", indx, list);

   RooArgSet normSet{x};

   indx.setIndex(0);
   EXPECT_EQ(pdf.getVal(normSet), gaus1.getVal(normSet));

   indx.setIndex(1);
   EXPECT_EQ(pdf.getVal(normSet), gaus2.getVal(normSet));
}

TEST(RooMultiPdfTest, FitConvergesAndReturnsReasonableResult)

{
   using namespace RooFit;
   const EvalBackend backend{EvalBackend::defaultValue()};

   RooRealVar x("x", "x", -10, 10);

   RooRealVar m1("mean1", "mean1", 0., -10, 10);
   RooRealVar s1("sigma1", "sigma1", 4., 0.001, 10.);
   RooRealVar m2("mean2", "mean2", 2., -10, 10);
   RooRealVar s2("sigma2", "sigma2", 4., 0.001, 10.);
   RooGaussian gaus1("gaus1", "gaus1", x, m1, s1);
   RooGaussian gaus2("gaus2", "gaus2", x, m2, s2);

   RooCategory indx("my_special_index", "my_index");
   indx.setConstant();

   RooMultiPdf multipdf("roomultipdf", "pdfs", indx, RooArgList{gaus1, gaus2});
   indx.setIndex(0); // choose first gaussian initially

   std::unique_ptr<RooDataSet> data{gaus1.generate(x, 1000)};

   // Fit 1 - RooMultiPdf fit

   std::unique_ptr<RooAbsReal> nll{multipdf.createNLL(*data, backend)};

   RooMinimizer minim{*nll};
   minim.setStrategy(0);
   int status = minim.minimize("Minuit2", "");

   // Resetting the initial parameters

   m1.setVal(0.);
   s1.setVal(4.);

   m1.setError(0.0);
   s1.setError(0.0);
   // Fit 2 - Reference fit
   std::unique_ptr<RooAbsReal> nll1{gaus1.createNLL(*data, backend)};

   RooMinimizer minim1{*nll1};
   minim1.setStrategy(0);
   int status1 = minim1.minimize("Minuit2", "");

   int n_param_gaus1 = countFloatingParametersIncludingObservable(gaus1);

   double first_fit = nll->getVal();
   double ref_fit = nll1->getVal() + 0.5 * n_param_gaus1; // 1.5 because the gaussian has 3 param*0.5
   // Now test the results
   EXPECT_EQ(status, 0) << "Fit 1 did not converge.";
   EXPECT_EQ(status1, 0) << "Fit 2 did not converge.";
   EXPECT_TRUE(std::isfinite(first_fit)) << "NLL is not finite.";
   EXPECT_TRUE(std::isfinite(ref_fit)) << "NLL1 is not finite.";
   EXPECT_DOUBLE_EQ(first_fit, ref_fit);

   // Check that the correct number of PDFs are present
   EXPECT_EQ(multipdf.getNumPdfs(), 2);

   //  check fitted parameter
   EXPECT_NEAR(m1.getVal(), 0.0, 0.2);
   EXPECT_NEAR(s1.getVal(), 4.0, 0.2);

   // Check whether RooMultiPdf chooses the correct index

   for (int i = 0; i < multipdf.getNumPdfs(); ++i) {

      indx.setIndex(i);

      std::unique_ptr<RooAbsReal> nll_multi{multipdf.createNLL(*data, backend)};

      RooAbsPdf *selectedPdf = multipdf.getPdf(i);
      std::unique_ptr<RooAbsReal> nll_direct{selectedPdf->createNLL(*data, backend)};

      int n_param = countFloatingParametersIncludingObservable(*selectedPdf);

      double multi = nll_multi->getVal();
      double direct = nll_direct->getVal() + 0.5 * n_param;

      std::cout << "PDF index " << i << ": n_param = " << n_param << ", direct+penalty = " << direct
                << ", multipdf = " << multi << std::endl;

      EXPECT_NEAR(multi, direct, 1e-6) << "Mismatch at index " << i;
   }
}

TEST(RooMultiPdfTest, PenaltyTermIsAppliedCorrectly)
{
   using namespace RooFit;
   const EvalBackend backend = EvalBackend::defaultValue();

   RooRealVar x("x", "x", -10, 10);

   RooRealVar mean("mean", "mean", 0, -5, 5);
   RooRealVar sigma("sigma", "sigma", 1, 0.1, 10);

   RooGaussian gauss1("gauss1", "gauss1", x, mean, sigma);
   RooGaussian gauss2("gauss2", "gauss2", x, mean, sigma);

   RooCategory index("index", "index");

   RooArgList pdfList(gauss1, gauss2);

   RooMultiPdf multiPdf("multiPdf", "multiPdf", index, pdfList);
   index.setConstant();

   index.setIndex(0);

   std::unique_ptr<RooDataSet> data{gauss1.generate(x, 100)};

   std::unique_ptr<RooAbsReal> nll_gauss1{gauss1.createNLL(*data, backend)};

   std::unique_ptr<RooAbsReal> nll_multi{multiPdf.createNLL(*data, backend)};

   double val_gauss1 = nll_gauss1->getVal();
   double val_multi = nll_multi->getVal();
   int n_params = countFloatingParametersIncludingObservable(gauss1);

   const double expected_penalty = 0.5 * n_params;
   const double delta = val_multi - val_gauss1;

   std::cout << "NLL(gauss1):     " << val_gauss1 << std::endl;
   std::cout << "NLL(multiPdf):   " << val_multi << std::endl;
   std::cout << "Expected penalty: " << expected_penalty << std::endl;
   std::cout << "Delta:           " << delta << std::endl;

   EXPECT_TRUE(std::isfinite(val_gauss1));
   EXPECT_TRUE(std::isfinite(val_multi));

   EXPECT_NEAR(delta, expected_penalty, 1e-6) << "Penalty term not correctly applied.";
}

// Test that the minimizer can correctly work even with disconnected floating
// parameters (it is expected to temporarily freeze them during the
// minimization).
TEST(RooMultiPdfTest, Minimization)
{
   RooRealVar x("x", "x", -10, 10);

   RooRealVar m1("mean1", "mean1", 0.);
   RooRealVar s1("sigma1", "sigma1", 4., 0.001, 10.);

   RooRealVar m2("mean2", "mean2", 0.5);
   RooRealVar s2("sigma2", "sigma2", 4., 0.001, 10.);

   RooGaussian gaus1("gaus1", "gaus1", x, m1, s1);
   RooGaussian gaus2("gaus2", "gaus2", x, m2, s2);

   RooCategory indx("my_special_index", "my_index");

   RooMultiPdf pdf("mult", "multi_pdf", indx, RooArgList{gaus1, gaus2});

   indx.setConstant();

   std::unique_ptr<RooAbsData> data{pdf.generate(x, 10000)};

   // Move parameter away from the value used to generate the dataset in order
   // to make the fit non-trivial.
   s1.setVal(3.);
   s2.setVal(3.);

   std::unique_ptr<RooAbsReal> nll1{gaus1.createNLL(*data)};
   std::unique_ptr<RooAbsReal> nll2{gaus2.createNLL(*data)};
   std::unique_ptr<RooAbsReal> nll{pdf.createNLL(*data)};

   const int nParams1 = 1 + 1; // plus one observable
   const int nParams2 = 1 + 1; // plus one observable

   int printLevel = -1;
   int nPdfs = 2;

   RooArgSet params{s1, s2};
   RooArgSet origParams;
   params.snapshot(origParams);

   RooMinimizer minim1{*nll1};
   minim1.setPrintLevel(printLevel);
   minim1.minimize("Minuit2", "");
   // Manually adding the penalty term
   double nllVal1 = nll1->getVal() + 0.5 * nParams1;
   params.assign(origParams);

   RooMinimizer minim2{*nll2};
   minim2.setPrintLevel(printLevel);
   minim2.minimize("Minuit2", "");
   double nllVal2 = nll2->getVal() + 0.5 * nParams2;
   params.assign(origParams);

   std::vector<double> multiNllVals;

   RooMinimizer minim{*nll};
   minim.setPrintLevel(printLevel);

   // Reuse the same minimizer to minimize for the different pdf choices one
   // after the other.
   for (int i = 0; i < nPdfs; ++i) {
      indx.setIndex(i);
      minim.minimize("Minuit2", "");
      multiNllVals.push_back(nll->getVal());
      params.assign(origParams);
   }

   // Validate the results
   EXPECT_DOUBLE_EQ(multiNllVals[0], nllVal1);
   EXPECT_DOUBLE_EQ(multiNllVals[1], nllVal2);
}
// When the index category is among the observables and the component pdfs are
// extended, generating from a RooMultiPdf should sample the index category
// according to the relative expected yields of the components, and generate the
// observables from the selected component (like RooSimultaneous does).
// Covers https://github.com/root-project/root/issues/22916
TEST(RooMultiPdfTest, ExtendedGenerationSamplesIndexByYields)
{
   using namespace RooFit;

   RooRandom::randomGenerator()->SetSeed(1337);

   RooRealVar x("x", "x", -10, 10);

   RooRealVar m1("m1", "m1", 0.0);
   RooRealVar s1("s1", "s1", 1.0);
   RooRealVar m2("m2", "m2", 3.0);
   RooRealVar s2("s2", "s2", 2.0);
   RooGaussian g1("g1", "g1", x, m1, s1);
   RooGaussian g2("g2", "g2", x, m2, s2);

   // Extend the components with different expected yields, so that the
   // category split is expected to be 100:300, i.e. 25% : 75%.
   RooRealVar n1("n1", "n1", 100.0);
   RooRealVar n2("n2", "n2", 300.0);
   RooExtendPdf e1("e1", "e1", g1, n1);
   RooExtendPdf e2("e2", "e2", g2, n2);

   RooCategory cat("cat", "cat");
   RooMultiPdf multiPdf("multiPdf", "multiPdf", cat, RooArgList{e1, e2});

   const std::string lab0 = cat.lookupName(0);
   const std::string lab1 = cat.lookupName(1);

   const int nEvents = 40000;
   std::unique_ptr<RooDataSet> data{multiPdf.generate({x, cat}, nEvents)};

   ASSERT_EQ(data->numEntries(), nEvents);

   const double n0 = data->sumEntries(("cat==cat::" + lab0).c_str());
   const double n1Obs = data->sumEntries(("cat==cat::" + lab1).c_str());

   ASSERT_EQ(n0 + n1Obs, nEvents);

   const double expectedFrac0 = 100. / 400.;
   const double expectedFrac1 = 300. / 400.;

   // Tolerance of 0.02 corresponds to roughly 9 sigma of the binomial
   // fluctuation for this sample size, so this is a robust check.
   EXPECT_NEAR(n0 / nEvents, expectedFrac0, 0.02) << "Index category not sampled proportionally to yields";
   EXPECT_NEAR(n1Obs / nEvents, expectedFrac1, 0.02) << "Index category not sampled proportionally to yields";

   // Verify that the observable values are actually generated from the
   // component associated with each category state.
   EXPECT_NEAR(data->mean(x, ("cat==cat::" + lab0).c_str()), m1.getVal(), 0.1);
   EXPECT_NEAR(data->mean(x, ("cat==cat::" + lab1).c_str()), m2.getVal(), 0.1);
}

// If the index category is not requested among the observables, the RooMultiPdf
// should generate from the currently-selected component pdf (unchanged legacy
// behaviour).
TEST(RooMultiPdfTest, GenerationWithoutIndexUsesCurrentPdf)
{
   RooRandom::randomGenerator()->SetSeed(1337);

   RooRealVar x("x", "x", -10, 10);

   RooRealVar m1("m1", "m1", -3.0);
   RooRealVar s1("s1", "s1", 0.5);
   RooRealVar m2("m2", "m2", 3.0);
   RooRealVar s2("s2", "s2", 0.5);
   RooGaussian g1("g1", "g1", x, m1, s1);
   RooGaussian g2("g2", "g2", x, m2, s2);

   RooCategory cat("cat", "cat");
   RooMultiPdf multiPdf("multiPdf", "multiPdf", cat, RooArgList{g1, g2});

   cat.setIndex(1); // select the second component

   std::unique_ptr<RooDataSet> data{multiPdf.generate(x, 10000)};

   // All events should come from the currently-selected component (g2).
   EXPECT_NEAR(data->mean(x), m2.getVal(), 0.1);
}

TEST(RooMultiReal, SelectsCorrectModel)
{
   RooRealVar x("x", "x", -10, 10);
   x.setVal(2.0);

   RooRealVar model1("model1", "model1", 5.0);
   RooRealVar model2("model2", "model2", 10.0);

   RooCategory indx("index", "index");

   RooArgList models{model1, model2};

   RooMultiReal multiReal("multiReal", "multi_real", indx, models);

   RooArgSet normSet{x};

   indx.setIndex(0);
   EXPECT_DOUBLE_EQ(multiReal.getVal(normSet), model1.getVal());

   indx.setIndex(1);
   EXPECT_DOUBLE_EQ(multiReal.getVal(normSet), model2.getVal());
}

TEST(RooMultiReal, EvaluateAndParameterAccess_Hook)
{
   RooRealVar x("x", "x", -10, 10);

   RooRealVar model1("model1", "model1", 5.0, 0., 10.);
   RooRealVar model2("model2", "model2", 10.0, 5., 15.);

   RooCategory indx("index", "index");

   RooMultiReal multiReal("multiReal", "multi_real", indx, RooArgList{model1, model2});

   indx.setIndex(0);
   EXPECT_EQ(multiReal.getVal(), model1.getVal());

   indx.setIndex(1);
   EXPECT_EQ(multiReal.getVal(), model2.getVal());

   // Prepare the observables
   RooArgSet observables(x);

   // Prepare an empty parameter list for getParametersHook to fill
   RooArgSet params;
   indx.setIndex(0);
   multiReal.getParameters(&observables, params, true);
   EXPECT_TRUE(params.find("model1") != nullptr);
   EXPECT_TRUE(params.find("model2") == nullptr);

   indx.setIndex(1);
   multiReal.getParameters(&observables, params, true);
   EXPECT_TRUE(params.find("model1") == nullptr);
   EXPECT_TRUE(params.find("model2") != nullptr);
}
TEST(RooMultiPdf, StripDisconnectedParameterTest)
{
   // Observable
   RooRealVar x("x", "x", -10, 10);

   // Parameters for two Gaussians
   RooRealVar mean1("mean1", "mean1", 0.0, -5.0, 5.0);
   RooRealVar sigma1("sigma1", "sigma1", 1.0, 0.1, 10.0);
   RooGaussian gauss1("gauss1", "gauss1", x, mean1, sigma1);

   RooRealVar mean2("mean2", "mean2", 2.0, -5.0, 5.0);
   RooRealVar sigma2("sigma2", "sigma2", 2.0, 0.1, 10.0);
   RooGaussian gauss2("gauss2", "gauss2", x, mean2, sigma2);

   // Category index
   RooCategory cat("cat", "cat");
   cat.defineType("first", 0);
   cat.defineType("second", 1);

   RooArgList pdfs(gauss1, gauss2);
   RooMultiPdf multiPdf("multiPdf", "multiPdf", cat, pdfs);

   RooArgSet observables(x);
   RooArgSet params;

   // --- Case 1: stripDisconnected = true ---
   cat.setIndex(0);
   multiPdf.getParameters(&observables, params, true);

   EXPECT_TRUE(params.find("mean1") != nullptr);
   EXPECT_TRUE(params.find("sigma1") != nullptr);

   EXPECT_TRUE(params.find("mean2") == nullptr);  // should be stripped out
   EXPECT_TRUE(params.find("sigma2") == nullptr); // should be stripped out

   // --- Case 2: stripDisconnected = false ---
   params.removeAll();
   cat.setIndex(0);
   multiPdf.getParameters(&observables, params, false);

   // Now both models' parameters should be present
   EXPECT_TRUE(params.find("mean1") != nullptr);
   EXPECT_TRUE(params.find("sigma1") != nullptr);
   EXPECT_TRUE(params.find("mean2") != nullptr);
   EXPECT_TRUE(params.find("sigma2") != nullptr);
}