## \file
## \ingroup tutorial_roofit_main
## \notebook -js
## Likelihood and minimization: approximate a likelihood with a simplified
## chi-square model for cheap toy studies.
##
## When the asymptotic approximation is not valid (for example in EFT
## analyses, where the parameters of interest enter the yields non-linearly),
## physics results have to be validated with toys. Running toys with the full
## model can be prohibitively expensive if it has many nuisance parameters. A
## common approach is to approximate the full likelihood with a much smaller
## chi-square model:
##
##   1. Fit the full model with one unconstrained yield parameter per bin
##      (like the ShapeFactors in HistFactory), so that the yields absorb the
##      parameters of interest.
##   2. Run Hesse to obtain the covariance matrix of the fitted yields, in
##      which the effect of all profiled nuisance parameters is absorbed.
##   3. Build a multivariate-Gaussian ("chi-square") model in which the bin
##      yields are expressed by an external parameterization with few degrees
##      of freedom, for example yields that are quadratic in the EFT
##      couplings.
##   4. Use the fast chi-square model for the toys.
##
## This tutorial demonstrates how steps 2. and 3. are automated by
## RooFitResult.createChi2Pdf() and RooFitResult.createChi2DataSet(), and
## validates the simplified model against a fit of the parameter of interest
## with the full model.
##
## \macro_image
## \macro_code
## \macro_output
##
## \date August 2026
## \author Jonas Rembser (CERN)

import ROOT

n_bins = 5

# Per-bin background yields and coefficients of the signal yields, which are
# quadratic in the EFT coupling "c":
#   yield_i(c) = bkg_i + s1_i * c + s2_i * c^2
bkg = [2000.0, 3000.0, 4000.0, 2500.0, 1500.0]
s1 = [50.0, 100.0, 200.0, 150.0, 100.0]
s2 = [10.0, 20.0, 40.0, 30.0, 20.0]

# The "full" model of the measurement: one Poisson counting experiment per
# bin. The yields are unconstrained parameters "mu_<bin>", shared among all
# EFT hypotheses, and they are multiplied by a common scale uncertainty
# implemented with the constrained nuisance parameter "alpha". In a real
# analysis, this model would have many more nuisance parameters, and
# evaluating it would be much more expensive than evaluating the simplified
# model derived from it below.
ws = ROOT.RooWorkspace("ws")
ws.factory("expr::kappa('1.0 + 0.05 * alpha', alpha[0, -5, 5])")
observables = ROOT.RooArgSet()
bin_params = ROOT.RooArgList()
for i in range(n_bins):
    ws.factory(f"Poisson::pois_{i}(n_{i}[0, 100000], prod::expYield_{i}(mu_{i}[{bkg[i]}, 1, 100000], kappa))")
    observables.add(ws[f"n_{i}"])
    bin_params.add(ws[f"mu_{i}"])
ws.factory("Gaussian::constraint(alpha, 0.0, 1.0)")
ws.factory("PROD::model(constraint," + ",".join(f"pois_{i}" for i in range(n_bins)) + ")")

# The external parameterization of the bin yields that will be used in the
# simplified model: quadratic functions of the coupling "c".
ws.factory("c[0, -10, 10]")
predictions = ROOT.RooArgList()
for i in range(n_bins):
    ws.factory(f"PolyVar::pred_{i}(c, {{{bkg[i]}, {s1[i]}, {s2[i]}}})")
    predictions.add(ws[f"pred_{i}"])

# Generate the "observed" dataset: one Poisson fluctuation of the yields
# expected for a true coupling of c = 2.
ROOT.RooRandom.randomGenerator().SetSeed(12345)
c_true = 2.0
for i in range(n_bins):
    ws[f"mu_{i}"].setVal(bkg[i] + s1[i] * c_true + s2[i] * c_true**2)
data = ws["model"].generate(observables, 1)

# Step 1: fit the full model with the unconstrained yield parameters. The
# covariance matrix from Hesse (run by default in fitTo()) absorbs the effect
# of the profiled nuisance parameter "alpha".
fit_result = ws["model"].fitTo(data, Save=True, PrintLevel=-1)
fit_result.Print()

# Steps 2. and 3.: build the simplified chi-square model from the fit result,
# replacing the measured yields by the quadratic parameterization, together
# with the single-entry dataset of the measured yields.
chi2_pdf = fit_result.createChi2Pdf(bin_params, predictions)
chi2_data = fit_result.createChi2DataSet(bin_params)

# Fit the coupling with the simplified model.
c = ws["c"]
fit_result_chi2 = chi2_pdf.fitTo(chi2_data, Save=True, PrintLevel=-1)
c_chi2 = c.getVal()
c_chi2_err = c.getError()

# Validation: fit the coupling directly with the full model, where the yields
# are given by the quadratic parameterization. The two results should agree
# well if the Gaussian approximation of the full likelihood is valid.
for i in range(n_bins):
    ws.factory(f"Poisson::poisFull_{i}(n_{i}, prod::expYieldFull_{i}(pred_{i}, kappa))")
ws.factory("PROD::modelFull(constraint," + ",".join(f"poisFull_{i}" for i in range(n_bins)) + ")")
# Seed the expensive full fit with the result of the cheap chi-square fit.
c.setVal(c_chi2)
c.setError(c_chi2_err)
ws["alpha"].setVal(0.0)
fit_result_full = ws["modelFull"].fitTo(data, Save=True, PrintLevel=-1)

print(f"\nFitted coupling with the full model       : {c.getVal():.4f} +/- {c.getError():.4f}")
print(f"Fitted coupling with the chi-square model : {c_chi2:.4f} +/- {c_chi2_err:.4f}\n")

# Step 4: run toys with the simplified model. Each toy is a single "event" in
# the pseudo-observables, i.e. one vector of measured yields. Here the toys
# are thrown at the best-fit coupling from the chi-square fit; the pull
# distribution of the coupling validates the parameter uncertainty. Small
# deviations from a unit Gaussian are expected, since the yields are
# non-linear in the coupling: quantifying such effects is exactly why toys are
# needed in the first place.
c.setVal(c_chi2)
c.setError(c_chi2_err)
chi2_observables = chi2_pdf.getObservables(chi2_data)
mc_study = ROOT.RooMCStudy(chi2_pdf, chi2_observables, ROOT.RooFit.Silence())
mc_study.generateAndFit(500, 1)

pull_frame = mc_study.plotPull(c, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(True))
pull_frame.SetTitle("Pull of the coupling in toys from the chi-square model")

canvas = ROOT.TCanvas("rf620_simplified_toys", "rf620_simplified_toys", 600, 600)
pull_frame.Draw()

canvas.SaveAs("rf620_simplified_toys.png")
