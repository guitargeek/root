/// \file
/// \ingroup tutorial_roostats
/// \notebook
/// Bayesian calculator: basic example
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Gregory Schott

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMsgService.h"

#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "TCanvas.h"

using namespace RooFit;
using namespace RooStats;

void rs701_BayesianCalculator(bool useBkg = true, double confLevel = 0.90)
{

   RooWorkspace *w = new RooWorkspace("w");
   w->factory("SUM::pdf(s[0.001,15]*Uniform(x[0,1]),b[1,0,2]*Uniform(x))");
   w->factory("Gaussian::prior_b(b,1,1)");
   w->factory("Uniform::priorPOI(s)");
   w->factory("GenericPdf::priorPOI2('1/sqrt(@0)',s)");

   w->factory("n[3]"); // observed number of events
   // create a data set with n observed events
   RooDataSet data("data", "", {*w->var("x"), *w->var("n")}, RooFit::WeightVar("n"));
   data.add({*(w->var("x"))}, w->var("n")->getVal());

   // to suppress messages when pdf goes to zero
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   // The model for a given prior of the parameter of interest s. The prior of
   // the nuisance parameter b is no part of the pdf: it multiplies the
   // likelihood only once, like a constraint term.
   auto makeModelConfig = [&](const char *name, const char *priorPOI) {
      auto modelConfig = std::make_unique<ModelConfig>(name, w);
      modelConfig->SetPdf("pdf");
      modelConfig->SetObservables(RooArgSet(*w->var("x")));
      modelConfig->SetParametersOfInterest(RooArgSet(*w->var("s")));
      if (useBkg)
         modelConfig->SetNuisanceParameters(RooArgSet(*w->var("b")));
      modelConfig->SetExternalConstraints(RooArgSet(*w->pdf("prior_b")));
      modelConfig->SetPriorPdf(priorPOI);
      return modelConfig;
   };

   double size = 1. - confLevel;
   std::cout << "\nBayesian Result using a Flat prior " << std::endl;
   std::unique_ptr<ModelConfig> modelConfig = makeModelConfig("modelConfig", "priorPOI");
   BayesianCalculator bcalc(data, *modelConfig);
   bcalc.SetTestSize(size);
   SimpleInterval *interval = bcalc.GetInterval();
   double cl = bcalc.ConfidenceLevel();
   std::cout << cl << "% CL central interval: [ " << interval->LowerLimit() << " - " << interval->UpperLimit()
             << " ] or " << cl + (1. - cl) / 2 << "% CL limits\n";
   RooPlot *plot = bcalc.GetPosteriorPlot();
   TCanvas *c1 = new TCanvas("c1", "Bayesian Calculator Result");
   c1->Divide(1, 2);
   c1->cd(1);
   plot->Draw();
   c1->Update();

   std::cout << "\nBayesian Result using a 1/sqrt(s) prior  " << std::endl;
   std::unique_ptr<ModelConfig> modelConfig2 = makeModelConfig("modelConfig2", "priorPOI2");
   BayesianCalculator bcalc2(data, *modelConfig2);
   bcalc2.SetTestSize(size);
   SimpleInterval *interval2 = bcalc2.GetInterval();
   cl = bcalc2.ConfidenceLevel();
   std::cout << cl << "% CL central interval: [ " << interval2->LowerLimit() << " - " << interval2->UpperLimit()
             << " ] or " << cl + (1. - cl) / 2 << "% CL limits\n";

   RooPlot *plot2 = bcalc2.GetPosteriorPlot();
   c1->cd(2);
   plot2->Draw();
   gPad->SetLogy();
   c1->Update();

   // observe one event while expecting one background event -> the 95% CL upper limit on s is 4.10
   // observe one event while expecting zero background event -> the 95% CL upper limit on s is 4.74
}
