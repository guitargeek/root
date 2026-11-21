## \file
## \ingroup tutorial_roofit_main
## \notebook -nodraw
## Likelihood and minimization: fitting with constraints
##
## \macro_code
## \macro_output
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Create model and dataset
# ----------------------------------------------

# Construct a Gaussian pdf
x = ROOT.RooRealVar("x", "x", -10, 10)

m = ROOT.RooRealVar("m", "m", 0, -10, 10)
s = ROOT.RooRealVar("s", "s", 2, 0.1, 10)
gauss = ROOT.RooGaussian("gauss", "gauss(x,m,s)", x, m, s)

# Construct a flat pdf (polynomial of 0th order)
poly = ROOT.RooPolynomial("poly", "poly(x)", x)

# model = f*gauss + (1-f)*poly
f = ROOT.RooRealVar("f", "f", 0.5, 0.0, 1.0)
model = ROOT.RooAddPdf("model", "model", [gauss, poly], [f])

# Generate small dataset for use in fitting below
d = model.generate(x, 50)

# Create constraint pdf
# -----------------------------------------

# Construct Gaussian constraint pdf on parameter f at 0.8 with
# resolution of 0.1
fconstraint = ROOT.RooGaussian("fconstraint", "fconstraint", f, 0.8, 0.1)

# A constraint term is no part of the pdf of the observable x: the pdf is the
# probability density for one event, while the constraint term multiplies the
# likelihood only once. Therefore, the constraint is not multiplied with the
# pdf in a RooProdPdf, but it is passed to the likelihood separately.

# Fit model (without use of constraint term)
r1 = model.fitTo(d, Save=True, PrintLevel=-1)

# Method 1 - specify constraint when fitting
# ------------------------------------------------------------------------------------------

# Fit model with constraint term on parameter f
r2 = model.fitTo(d, ExternalConstraints=[fconstraint], Save=True, PrintLevel=-1)

# Method 2 - keep model and constraint in a ModelConfig
# ------------------------------------------------------------------------------------------

# The recommended way to keep a pdf together with its constraint terms is a
# ModelConfig in a RooWorkspace: you can't forget to pass the constraints to
# the fit, and the RooStats tools use them to compute limits and intervals.

# Construct another Gaussian constraint pdf on parameter f at 0.2 with
# resolution of 0.1
fconstext = ROOT.RooGaussian("fconstext", "fconstext", f, 0.2, 0.1)

ws = ROOT.RooWorkspace("ws")
modelConfig = ROOT.RooFit.ModelConfig("ModelConfig", ws)
# The ModelConfig imports the objects into the workspace
modelConfig.SetPdf(model)
modelConfig.SetExternalConstraints([fconstext])
modelConfig.SetObservables("x")
modelConfig.SetParametersOfInterest("f")

# ModelConfig.fitTo() passes the constraint terms to the fit
r3 = modelConfig.fitTo(d, Save=True, PrintLevel=-1)

# Print the fit results
print("fit result without constraint (data generated at f=0.5)")
r1.Print("v")
print("fit result with constraint (data generated at f=0.5, is f=0.8+/-0.1)")
r2.Print("v")
print("fit result with (another) constraint in a ModelConfig (data generated at f=0.5, is f=0.2+/-0.1)")
r3.Print("v")
