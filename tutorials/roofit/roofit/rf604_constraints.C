/// \file
/// \ingroup tutorial_roofit_main
/// \notebook -nodraw
/// Likelihood and minimization: fitting with constraints
///
/// \macro_code
/// \macro_output
///
/// \date 07/2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooFit/ModelConfig.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf604_constraints()
{

   // C r e a t e   m o d e l  a n d   d a t a s e t
   // ----------------------------------------------

   // Construct a Gaussian pdf
   RooRealVar x("x", "x", -10, 10);

   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 2, 0.1, 10);
   RooGaussian gauss("gauss", "gauss(x,m,s)", x, m, s);

   // Construct a flat pdf (polynomial of 0th order)
   RooPolynomial poly("poly", "poly(x)", x);

   // Construct model = f*gauss + (1-f)*poly
   RooRealVar f("f", "f", 0.5, 0., 1.);
   RooAddPdf model("model", "model", RooArgSet(gauss, poly), f);

   // Generate small dataset for use in fitting below
   std::unique_ptr<RooDataSet> d{model.generate(x, 50)};

   // C r e a t e   c o n s t r a i n t   p d f
   // -----------------------------------------

   // Construct Gaussian constraint pdf on parameter f at 0.8 with resolution of 0.1
   RooGaussian fconstraint("fconstraint", "fconstraint", f, 0.8, 0.2);

   // A constraint term is no part of the pdf of the observable x: the pdf is the
   // probability density for one event, while the constraint term multiplies the
   // likelihood only once. Therefore, the constraint is not multiplied with the
   // pdf in a RooProdPdf, but it is passed to the likelihood separately.

   // Fit model (without use of constraint term)
   std::unique_ptr<RooFitResult> r1{model.fitTo(*d, Save(), PrintLevel(-1))};

   // M E T H O D   1   -   S p e c i f y   c o n s t r a i n t   w h e n   f i t t i n g
   // -----------------------------------------------------------------------------------

   // Fit model with constraint term on parameter f
   std::unique_ptr<RooFitResult> r2{model.fitTo(*d, ExternalConstraints(fconstraint), Save(), PrintLevel(-1))};

   // M E T H O D   2   -   K e e p   m o d e l   a n d   c o n s t r a i n t   i n   a   M o d e l C o n f i g
   // -----------------------------------------------------------------------------------------------------------

   // The recommended way to keep a pdf together with its constraint terms is a
   // ModelConfig in a RooWorkspace: you can't forget to pass the constraints to
   // the fit, and the RooStats tools use them to compute limits and intervals.

   // Construct another Gaussian constraint pdf on parameter f at 0.2 with resolution of 0.1
   RooGaussian fconstext("fconstext", "fconstext", f, 0.2, 0.1);

   RooWorkspace ws("ws");
   RooFit::ModelConfig modelConfig("ModelConfig", &ws);
   // The ModelConfig imports the objects into the workspace
   modelConfig.SetPdf(model);
   modelConfig.SetExternalConstraints(fconstext);
   modelConfig.SetObservables("x");
   modelConfig.SetParametersOfInterest("f");

   // ModelConfig::fitTo() passes the constraint terms to the fit
   std::unique_ptr<RooFitResult> r3{modelConfig.fitTo(*d, Save(), PrintLevel(-1))};

   // Print the fit results
   cout << "fit result without constraint (data generated at f=0.5)" << endl;
   r1->Print("v");
   cout << "fit result with constraint (data generated at f=0.5, constraint is f=0.8+/-0.2)" << endl;
   r2->Print("v");
   cout << "fit result with (another) constraint in a ModelConfig (data generated at f=0.5, constraint is f=0.2+/-0.1)"
        << endl;
   r3->Print("v");
}
