/// \file
/// \ingroup tutorial_roofit_main
/// \notebook -js
/// Likelihood and minimization: approximate a likelihood with a simplified
/// chi-square model for cheap toy studies.
///
/// When the asymptotic approximation is not valid (for example in EFT
/// analyses, where the parameters of interest enter the yields non-linearly),
/// physics results have to be validated with toys. Running toys with the full
/// model can be prohibitively expensive if it has many nuisance parameters. A
/// common approach is to approximate the full likelihood with a much smaller
/// chi-square model:
///
///   1. Fit the full model with one unconstrained yield parameter per bin
///      (like the ShapeFactors in HistFactory), so that the yields absorb the
///      parameters of interest.
///   2. Run Hesse to obtain the covariance matrix of the fitted yields, in
///      which the effect of all profiled nuisance parameters is absorbed.
///   3. Build a multivariate-Gaussian ("chi-square") model in which the bin
///      yields are expressed by an external parameterization with few degrees
///      of freedom, for example yields that are quadratic in the EFT
///      couplings.
///   4. Use the fast chi-square model for the toys.
///
/// This tutorial demonstrates how steps 2. and 3. are automated by
/// RooFitResult::createChi2Pdf() and RooFitResult::createChi2DataSet(), and
/// validates the simplified model against a fit of the parameter of interest
/// with the full model.
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \date August 2026
/// \author Jonas Rembser (CERN)

#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMCStudy.h>
#include <RooPlot.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <TCanvas.h>

#include <memory>

void rf620_simplified_toys()
{
   using namespace RooFit;

   const int nBins = 5;

   // Per-bin background yields and coefficients of the signal yields, which
   // are quadratic in the EFT coupling "c":
   //   yield_i(c) = bkg_i + s1_i * c + s2_i * c^2
   const double bkg[nBins] = {2000., 3000., 4000., 2500., 1500.};
   const double s1[nBins] = {50., 100., 200., 150., 100.};
   const double s2[nBins] = {10., 20., 40., 30., 20.};

   // The "full" model of the measurement: one Poisson counting experiment per
   // bin. The yields are unconstrained parameters "mu_<bin>", shared among
   // all EFT hypotheses, and they are multiplied by a common scale
   // uncertainty implemented with the constrained nuisance parameter "alpha".
   // In a real analysis, this model would have many more nuisance parameters,
   // and evaluating it would be much more expensive than evaluating the
   // simplified model derived from it below.
   RooWorkspace ws{"ws"};
   ws.factory("expr::kappa('1.0 + 0.05 * alpha', alpha[0, -5, 5])");
   RooArgSet observables;
   RooArgList binParams;
   std::string prodFactoryString = "PROD::model(constraint";
   for (int i = 0; i < nBins; ++i) {
      ws.factory(
         Form("Poisson::pois_%d(n_%d[0, 100000], prod::expYield_%d(mu_%d[%g, 1, 100000], kappa))", i, i, i, i, bkg[i]));
      observables.add(*ws.var(Form("n_%d", i)));
      binParams.add(*ws.var(Form("mu_%d", i)));
      prodFactoryString += Form(",pois_%d", i);
   }
   ws.factory("Gaussian::constraint(alpha, 0.0, 1.0)");
   ws.factory((prodFactoryString + ")").c_str());

   // The external parameterization of the bin yields that will be used in the
   // simplified model: quadratic functions of the coupling "c".
   ws.factory("c[0, -10, 10]");
   RooArgList predictions;
   for (int i = 0; i < nBins; ++i) {
      ws.factory(Form("PolyVar::pred_%d(c, {%g, %g, %g})", i, bkg[i], s1[i], s2[i]));
      predictions.add(*ws.function(Form("pred_%d", i)));
   }

   // Generate the "observed" dataset: one Poisson fluctuation of the yields
   // expected for a true coupling of c = 2.
   RooRandom::randomGenerator()->SetSeed(12345);
   const double cTrue = 2.;
   for (int i = 0; i < nBins; ++i) {
      ws.var(Form("mu_%d", i))->setVal(bkg[i] + s1[i] * cTrue + s2[i] * cTrue * cTrue);
   }
   std::unique_ptr<RooDataSet> data{ws.pdf("model")->generate(observables, 1)};

   // Step 1: fit the full model with the unconstrained yield parameters. The
   // covariance matrix from Hesse (run by default in fitTo()) absorbs the
   // effect of the profiled nuisance parameter "alpha".
   std::unique_ptr<RooFitResult> fr{ws.pdf("model")->fitTo(*data, Save(), PrintLevel(-1))};
   fr->Print();

   // Steps 2. and 3.: build the simplified chi-square model from the fit
   // result, replacing the measured yields by the quadratic parameterization,
   // together with the single-entry dataset of the measured yields.
   std::unique_ptr<RooAbsPdf> chi2Pdf = fr->createChi2Pdf(binParams, predictions);
   std::unique_ptr<RooDataSet> chi2Data = fr->createChi2DataSet(binParams);

   // Fit the coupling with the simplified model.
   RooRealVar &c = *ws.var("c");
   std::unique_ptr<RooFitResult> frChi2{chi2Pdf->fitTo(*chi2Data, Save(), PrintLevel(-1))};
   const double cChi2 = c.getVal();
   const double cChi2Err = c.getError();

   // Validation: fit the coupling directly with the full model, where the
   // yields are given by the quadratic parameterization. The two results
   // should agree well if the Gaussian approximation of the full likelihood
   // is valid.
   std::string prodFullFactoryString = "PROD::modelFull(constraint";
   for (int i = 0; i < nBins; ++i) {
      ws.factory(Form("Poisson::poisFull_%d(n_%d, prod::expYieldFull_%d(pred_%d, kappa))", i, i, i, i));
      prodFullFactoryString += Form(",poisFull_%d", i);
   }
   ws.factory((prodFullFactoryString + ")").c_str());
   // Seed the expensive full fit with the result of the cheap chi-square fit.
   c.setVal(cChi2);
   c.setError(cChi2Err);
   ws.var("alpha")->setVal(0.);
   std::unique_ptr<RooFitResult> frFull{ws.pdf("modelFull")->fitTo(*data, Save(), PrintLevel(-1))};

   std::cout << "\nFitted coupling with the full model       : " << c.getVal() << " +/- " << c.getError() << "\n";
   std::cout << "Fitted coupling with the chi-square model : " << cChi2 << " +/- " << cChi2Err << "\n\n";

   // Step 4: run toys with the simplified model. Each toy is a single "event"
   // in the pseudo-observables, i.e. one vector of measured yields. Here the
   // toys are thrown at the best-fit coupling from the chi-square fit; the
   // pull distribution of the coupling validates the parameter uncertainty.
   // Small deviations from a unit Gaussian are expected, since the yields are
   // non-linear in the coupling: quantifying such effects is exactly why toys
   // are needed in the first place.
   c.setVal(cChi2);
   c.setError(cChi2Err);
   std::unique_ptr<RooArgSet> chi2Observables{chi2Pdf->getObservables(*chi2Data)};
   RooMCStudy mcStudy{*chi2Pdf, *chi2Observables, Silence()};
   mcStudy.generateAndFit(500, 1);

   RooPlot *pullFrame = mcStudy.plotPull(c, Bins(20), FitGauss(true));
   pullFrame->SetTitle("Pull of the coupling in toys from the chi-square model");

   auto canvas = new TCanvas("rf620_simplified_toys", "rf620_simplified_toys", 600, 600);
   pullFrame->Draw();
}
