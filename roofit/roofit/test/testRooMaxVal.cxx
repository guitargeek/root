// Regression tests for RooAbsReal::getMaxVal / maxVal
// Covers the classes flagged in ROOT issue #12317. The contract is:
// when getMaxVal(nset) returns a non-zero code, maxVal(code) must be a
// valid upper bound on getVal(nset) for PDFs (on getVal() for plain
// RooAbsReal).

#include <RooBinSamplingPdf.h>
#include <RooCBShape.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistFunc.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooRealVar.h>
#include <RooWrapperPdf.h>

#include "gtest/gtest.h"

namespace {

// If getMaxVal returns a non-zero code, assert that maxVal is a valid
// upper bound on getVal(nset). If the code is zero, just return.
void expectUpperBound(RooAbsReal const &func, RooArgSet const &nset)
{
   const int code = func.getMaxVal(nset);
   if (code == 0) {
      return;
   }
   const double peak = func.getVal(nset);
   const double maxV = func.maxVal(code);
   EXPECT_GE(maxV, peak) << func.ClassName() << " maxVal(" << maxV << ") must be >= getVal(" << peak
                         << ") for nset of size " << nset.size();
}

} // namespace

TEST(RooMaxVal, UpperBoundContract)
{
   // Matches the reproducer in issue #12317.
   RooRealVar x{"x", "x", 2.5, 0, 5};
   RooRealVar y{"y", "y", 2.5, 0, 5};
   x.setBins(25);
   y.setBins(25);

   RooRealVar x0("x0", "x0", 2.5, -5., 5.);
   RooRealVar sigma("sigma", "sigma", 0.02, 1.E-4, 1.);
   RooRealVar alpha("a", "a", 1., 1.E-6, 100.);
   RooRealVar n("n", "n", 1., 1.E-6, 100.);
   RooCrystalBall crystalBall("cb1", "cb1", x, x0, sigma, alpha, n);
   RooCBShape cbShape("cb2", "cb2", x, x0, sigma, alpha, n);

   RooBinSamplingPdf binSamplingPdf("binSamplingPdf", "binSamplingPdf", x, crystalBall);

   RooDataHist templHist{"templHist", "templHist", x};
   templHist.set(x.getBin(), 100.0, -1.0);

   RooDataSet templData{"templData", "templData", x};
   for (int i = 0; i < x.numBins(); ++i) {
      x.setBin(i);
      templData.add(x);
   }
   x.setVal(2.5);

   RooHistFunc histFunc{"histFunc", "histFunc", x, templHist, 0};
   RooHistPdf histPdf{"histPdf", "histPdf", x, templHist, 0};
   RooWrapperPdf wrapperPdf{"wrapperPdf", "wrapperPdf", histFunc};
   RooKeysPdf keysPdf{"keysPdf", "keysPdf", x, templData};

   const RooArgSet normSetFull{x};
   const RooArgSet normSetExtra{x, y};

   // Prime any caches for the wrapper.
   wrapperPdf.getVal(normSetFull);

   for (auto const *nset : {&normSetFull, &normSetExtra}) {
      expectUpperBound(histFunc, *nset);
      expectUpperBound(histPdf, *nset);
      expectUpperBound(keysPdf, *nset);
      expectUpperBound(crystalBall, *nset);
      expectUpperBound(cbShape, *nset);
      expectUpperBound(wrapperPdf, *nset);
      expectUpperBound(binSamplingPdf, *nset);
   }
}

// Specifically check the bug fixed for RooHistPdf: maxVal used to return
// the raw bin weight * 1.05 (here 105), but the normalized PDF peaks at
// 5, so maxVal must now be close to 5 rather than 105.
TEST(RooMaxVal, HistPdfReturnsNormalizedMax)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   x.setBins(25);

   RooDataHist templHist{"templHist", "templHist", x};
   templHist.set(x.getBin(), 100.0, -1.0);
   x.setVal(2.5);

   RooHistPdf histPdf{"histPdf", "histPdf", x, templHist, 0};

   const RooArgSet nset{x};
   const int code = histPdf.getMaxVal(nset);
   ASSERT_NE(code, 0);

   const double peakVal = histPdf.getVal(nset);
   const double maxV = histPdf.maxVal(code);
   EXPECT_GE(maxV, peakVal);
   // Tight bound: within 10% of the true peak.
   EXPECT_LT(maxV, peakVal * 1.10);
}
