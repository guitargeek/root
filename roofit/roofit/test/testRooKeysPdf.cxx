// Tests for the RooKeysPdf and friends
// Authors: Jonas Rembser, CERN  07/2022

#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooKeysPdf.h>
#include <RooNDKeysPdf.h>
#include <RooNumIntConfig.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include "gtest/gtest.h"

// Test the support of RooKeysPdf and RooNDKeysPdf for weighted datasets.
TEST(RooKeysPdf, WeightedDataset)
{
   // We create data with 100 events at x = 5 and 400 events at x = 15. One
   // version will have 500 unweighted entries, the other will have only 2
   // entries with the weights 100 and 400 to represent the same data. The
   // resulting RooKeysPdfs should be identical for both datasets. Checking
   // this validates that dataset weights are correctly dealt with.

   RooRealVar x("x", "x", 0, 20);

   const std::size_t nEvents0 = 100;
   const std::size_t nEvents1 = 400;

   RooDataSet data1{"data1", "data1", x};
   RooDataSet data2{"data2", "data2", x, RooFit::WeightVar()};

   x.setVal(5);
   for (std::size_t i = 0; i < nEvents0; ++i) {
      data1.add(x, 1.0);
   }
   data2.add(x, double(nEvents0));
   x.setVal(15);
   for (std::size_t i = 0; i < nEvents1; ++i) {
      data1.add(x, 1.0);
   }
   data2.add(x, double(nEvents1));

   // Creating RooKeysPdf and RooNDKeysPdf with adaptive kernel and no
   // mirroring for both weighted and unweighted datasets.
   RooKeysPdf kest1("kest1", "kest1", x, data1, RooKeysPdf::NoMirror);
   RooKeysPdf kest2("kest2", "kest2", x, data2, RooKeysPdf::NoMirror);
   RooNDKeysPdf kestND1("kestND1", "kestND1", x, data1, "a");
   RooNDKeysPdf kestND2("kestND2", "kestND2", x, data2, "a");

   RooArgSet normSet{x};

   // Check if values for the unweighted and weighted datasets are the same
   double xVal = x.getMin();
   while (xVal < x.getMax()) {
      EXPECT_FLOAT_EQ(kest1.getVal(normSet), kest2.getVal(normSet));
      EXPECT_FLOAT_EQ(kestND1.getVal(normSet), kestND2.getVal(normSet));
      x.setVal(xVal);
      xVal += 0.1;
   }
}

// Test generation with proto data, covering GitHub issue #12286.
TEST(RooKeysPdf, GenerationWithProtoData)
{
   using namespace RooFit;

   RooRealVar x{"x", "", 0, 1};
   RooGenericPdf pdfX{"pdf_x", "x", {x}};

   std::unique_ptr<RooDataSet> dtBase{pdfX.generate(x, 10000)};

   RooKeysPdf pdfKeys{"pdf_keys", "", x, *dtBase, RooKeysPdf::MirrorBoth};

   RooRealVar y{"y", "", 0, 1};
   RooDataSet proto{"proto_y", "", y};
   proto.add(y);

   std::unique_ptr<RooDataSet> dtKeysWithProto{pdfKeys.generate(x, NumEvents(10000), ProtoData(proto))};

   std::unique_ptr<RooPlot> frame{x.frame()};
   dtKeysWithProto->plotOn(frame.get());
   pdfKeys.plotOn(frame.get());

   // If the dataset generation worked, the chi-square is not too terrible
   EXPECT_LE(frame->chiSquare(), 2.0);
}

// Test the analytic integral of RooNDKeysPdf over a named range, covering
// GitHub issue #6557. Analytic integration over a range is only exact when
// the kernels are not rotated into a decorrelated eigenbasis, so it is only
// enabled in that case (1D, or explicitly disabling rotation).
TEST(RooNDKeysPdf, RangedAnalyticIntegral1D)
{
   RooRealVar x("x", "x", -5, 15);

   RooDataSet data("data", "data", x);
   for (double xVal : {1.0, 2.5, 3.3, 4.8, 6.6, 7.8, 9.1, 10.0, 11.4, 13.2}) {
      x.setVal(xVal);
      data.add(x);
   }

   RooNDKeysPdf pdf("pdf", "pdf", RooArgList{x}, data, "a");

   x.setRange("sub", -2, 6);

   RooNDKeysPdf pdfNumInt(pdf);
   pdfNumInt.forceNumInt(true);

   std::unique_ptr<RooAbsReal> integral{pdf.createIntegral(x, "sub")};
   std::unique_ptr<RooAbsReal> numInt{pdfNumInt.createIntegral(x, "sub")};

   // getAnalyticalIntegral should have handed out the analytic code, i.e.
   // `integral` should not just be falling back to numeric integration.
   RooArgSet allVars{x};
   RooArgSet analVars;
   EXPECT_EQ(pdf.getAnalyticalIntegral(allVars, analVars, "sub"), 1);

   EXPECT_NEAR(integral->getVal(), numInt->getVal(), 1e-3 * numInt->getVal());
}

TEST(RooNDKeysPdf, RangedAnalyticIntegral2DUnrotated)
{
   RooRealVar x("x", "x", -5, 15);
   RooRealVar y("y", "y", -5, 15);
   RooArgSet vars{x, y};

   RooDataSet data("data", "data", vars);
   for (auto const &p : {std::pair{1.0, 2.0}, std::pair{2.5, 1.5}, std::pair{3.3, 4.1}, std::pair{4.8, 3.0},
                          std::pair{6.6, 5.5}, std::pair{7.8, 6.9}, std::pair{9.1, 8.0}, std::pair{10.0, 9.4},
                          std::pair{11.4, 10.1}, std::pair{13.2, 12.0}}) {
      x.setVal(p.first);
      y.setVal(p.second);
      data.add(vars);
   }

   // rotate = false: the kernels stay axis-aligned with x and y, so the
   // analytic ranged integral is exact.
   RooNDKeysPdf pdf("pdf", "pdf", RooArgList{x, y}, data, "a", 1.0, 3, /*rotate=*/false);

   x.setRange("sub", -2, 6);
   y.setRange("sub", -3, 5);

   RooNDKeysPdf pdfNumInt(pdf);
   pdfNumInt.forceNumInt(true);

   std::unique_ptr<RooAbsReal> integral{pdf.createIntegral(vars, "sub")};
   std::unique_ptr<RooAbsReal> numInt{pdfNumInt.createIntegral(vars, "sub")};

   EXPECT_NEAR(integral->getVal(), numInt->getVal(), 1e-2 * numInt->getVal());
}

// With rotation enabled in more than 1D, there is no simple closed form for
// the box integral, so we must keep falling back to numeric integration.
TEST(RooNDKeysPdf, RangedAnalyticIntegral2DRotatedFallsBackToNumeric)
{
   RooRealVar x("x", "x", -5, 15);
   RooRealVar y("y", "y", -5, 15);
   RooArgSet vars{x, y};

   RooDataSet data("data", "data", vars);
   for (auto const &p : {std::pair{1.0, 2.0}, std::pair{2.5, 1.5}, std::pair{3.3, 4.1}}) {
      x.setVal(p.first);
      y.setVal(p.second);
      data.add(vars);
   }

   RooNDKeysPdf pdf("pdf", "pdf", RooArgList{x, y}, data, "a", 1.0, 3, /*rotate=*/true);

   x.setRange("sub", -2, 6);
   y.setRange("sub", -3, 5);

   RooArgSet allVars{x, y};
   RooArgSet analVars;
   EXPECT_EQ(pdf.getAnalyticalIntegral(allVars, analVars, "sub"), 0);
}
