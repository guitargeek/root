/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2024
 *
 * Copyright (c) 2024, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <RooCmdArg.h>
#include <RooCurve.h>
#include <RooDataHist.h>
#include <RooGenericPdf.h>
#include <RooHelpers.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <Math/Util.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "gtest_wrapper.h"

// Cross-check to make sure the integration works correctly even if there is
// only one midpoint on the RooCurve. Covers GitHub issue #9838 (the reproducer
// in that issue was translated to this test).
TEST(RooPlot, Average)
{
   // Silence the info about numeric integration because we don't care about it
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING, 0u, RooFit::NumericIntegration, true};

   RooRealVar x("x", "x", 0, 50);
   RooGenericPdf func("func", "Test Function", "x", x);

   std::unique_ptr<RooPlot> xframe{x.frame()};

   func.plotOn(xframe.get(), RooFit::Name("funcCurve"));

   RooCurve *funcCurve = xframe->getCurve("funcCurve");

   const double tol = 1e-10;

   for (double i = 10; i < 11; i += 0.1) {
      double avg = funcCurve->average(i, i + 0.1);

      double xFirst = funcCurve->interpolate(i, tol);
      double xLast = funcCurve->interpolate(i + 0.1, tol);

      EXPECT_NEAR(avg, 0.5 * (xLast + xFirst), tol);
   }
}

/// Make sure that RooPlot::chiSquare() doesn't silently return NaN when the
/// data was plotted with a data error type that leaves the data points without
/// errors, like RooAbsData::Expected or RooAbsData::None. In that case, the
/// chi-square is not defined (and the bins must not be skipped, because that
/// would bias the test statistic), so an error is reported and NaN is returned
/// on purpose. Covers GitHub issue #21697.
TEST(RooPlot, ChiSquareWithoutDataErrors)
{
   using namespace RooFit;

   // Silence the info about numeric integration because we don't care about it
   RooHelpers::LocalChangeMsgLevel chmsglvl{RooFit::WARNING, 0u, RooFit::NumericIntegration, true};

   // Note the empty bins in the tails: they are skipped by the chi-square, so
   // they must not be counted as bins without errors either.
   const std::vector<double> binContent{0, 0, 0, 0, 0, 1, 3, 16, 17, 14, 18, 18, 7, 5, 1, 0, 0, 0, 0, 0};
   const int nBins = static_cast<int>(binContent.size());
   const int nNonEmptyBins =
      static_cast<int>(std::count_if(binContent.begin(), binContent.end(), [](double y) { return y != 0.0; }));

   TH1D h1{"h1", "h1", nBins, 0.0, 1.0};
   for (int i = 0; i < nBins; ++i) {
      h1.SetBinContent(i + 1, binContent[i]);
   }

   RooRealVar x{"x", "x", 0.0, 1.0};
   x.setBins(nBins);
   RooDataHist dataHist{"dataHist", "dataHist", RooArgList{x}, &h1};

   // Roughly the best-fit Gaussian for the data above (RooGenericPdf is used to
   // avoid a dependency on the RooFit library in this test).
   RooRealVar mean{"mean", "mean", 0.53};
   RooRealVar sigma{"sigma", "sigma", 0.11};
   RooGenericPdf pdf{"pdf", "pdf", "std::exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma))", {x, mean, sigma}};

   auto chiSquare = [&](RooAbsData::ErrorType etype, std::string &plotOnErrors, std::string &chi2Errors) {
      std::unique_ptr<RooPlot> frame{x.frame()};
      {
         RooHelpers::HijackMessageStream hijack{RooFit::ERROR, RooFit::InputArguments};
         dataHist.plotOn(frame.get(), DataError(etype), Name("data"));
         plotOnErrors = hijack.str();
      }
      pdf.plotOn(frame.get(), Name("curve"));
      RooHelpers::HijackMessageStream hijack{RooFit::ERROR, RooFit::Plotting};
      double out = frame->chiSquare("curve", "data", 2);
      chi2Errors = hijack.str();
      return out;
   };

   std::string plotOnErrors;
   std::string chi2Errors;

   // The supported error types still give a finite chi-square, without any error message
   for (auto etype : {RooAbsData::Poisson, RooAbsData::SumW2, RooAbsData::Auto}) {
      const double chi2 = chiSquare(etype, plotOnErrors, chi2Errors);
      EXPECT_TRUE(std::isfinite(chi2)) << "for error type " << etype;
      EXPECT_GT(chi2, 0.0) << "for error type " << etype;
      EXPECT_TRUE(plotOnErrors.empty()) << "for error type " << etype << ": " << plotOnErrors;
      EXPECT_TRUE(chi2Errors.empty()) << "for error type " << etype << ": " << chi2Errors;
   }

   // Only the non-empty bins may be reported as bins without errors: the empty
   // ones don't contribute to the chi-square in the first place.
   const std::string expectedCount = "has " + std::to_string(nNonEmptyBins) + " non-empty bin(s) with a zero error";

   // Expected errors can't be attached to the data points, which should already
   // be flagged by RooAbsData::plotOn()
   {
      const double chi2 = chiSquare(RooAbsData::Expected, plotOnErrors, chi2Errors);
      EXPECT_TRUE(std::isnan(chi2));
      EXPECT_NE(plotOnErrors.find("DataError(RooAbsData::Expected) is not supported"), std::string::npos)
         << plotOnErrors;
      EXPECT_NE(chi2Errors.find("chi-square is not defined"), std::string::npos) << chi2Errors;
      EXPECT_NE(chi2Errors.find(expectedCount), std::string::npos) << chi2Errors;
   }

   // Same for data plotted without any errors at all
   {
      const double chi2 = chiSquare(RooAbsData::None, plotOnErrors, chi2Errors);
      EXPECT_TRUE(std::isnan(chi2));
      EXPECT_TRUE(plotOnErrors.empty()) << plotOnErrors;
      EXPECT_NE(chi2Errors.find("chi-square is not defined"), std::string::npos) << chi2Errors;
      EXPECT_NE(chi2Errors.find(expectedCount), std::string::npos) << chi2Errors;
   }

   // Finally, pin down the chi-square definition that is documented in
   // RooCurve::chiSquare(): the numerator is the average of the plotted curve
   // over the bin, the denominator is the asymmetric error of the data point
   // (the lower one if the data is above the curve, the upper one otherwise),
   // empty and out-of-range bins are skipped and don't count towards the ndf,
   // and the returned value is chi2/ndf.
   {
      std::unique_ptr<RooPlot> frame{x.frame()};
      dataHist.plotOn(frame.get(), DataError(RooAbsData::Poisson), Name("data"));
      pdf.plotOn(frame.get(), Name("curve"));

      RooHist &hist = *frame->getHist("data");
      RooCurve &curve = *frame->getCurve("curve");

      const double xStart = curve.GetPointX(0);
      const double xStop = curve.GetPointX(curve.GetN() - 1);

      ROOT::Math::KahanSum<double> chi2Ref;
      int nRefBins = 0;
      for (int i = 0; i < hist.GetN(); ++i) {
         const double xi = hist.GetPointX(i);
         const double yi = hist.GetPointY(i);
         if (xi < xStart || xi > xStop || yi == 0.0) {
            continue;
         }
         const double avg = curve.average(xi - hist.GetEXlow()[i], xi + hist.GetEXhigh()[i]);
         const double err = yi > avg ? hist.GetEYlow()[i] : hist.GetEYhigh()[i];
         EXPECT_GT(err, 0.0) << "for bin " << i;
         const double pull = (yi - avg) / err;
         chi2Ref += pull * pull;
         ++nRefBins;
      }

      EXPECT_EQ(nRefBins, nNonEmptyBins);
      EXPECT_DOUBLE_EQ(frame->chiSquare("curve", "data", 2), chi2Ref.Sum() / (nRefBins - 2));
   }
}
