// Tests for the RooViewDataStore
// Authors: Jonas Rembser, CERN  07/2021

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooViewDataStore.h"
#include "RooMsgService.h"

#include "TRandom3.h"

#include "gtest/gtest.h"

#include <memory>

namespace {

void setupRooMsgService()
{
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
   RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
}

} // namespace

TEST(RooViewDataStore, UnweightedDataSetOneDim)
{
   const double xMin = -10.0;
   const double xMax = 20.0;
   const double meanVal = 5.0;
   const double sigmaVal = 2.0;
   const double meanStartVal = 2.0;
   const double sigmaStartVal = 3.0;
   const int randomSeed = 65539;
   const std::size_t nEntries = 100;

   using namespace RooFit;

   setupRooMsgService();

   TRandom3 rndm{randomSeed};

   std::vector<double> xValues;
   xValues.resize(nEntries);

   for (std::size_t i = 0; i < nEntries; ++i) {
      xValues[i] = rndm.Gaus(meanVal, sigmaVal);
      xValues[i] = std::max(xMin, xValues[i]);
      xValues[i] = std::min(xMax, xValues[i]);
   }

   RooRealVar x("x", "x", (xMax - xMin) / 2., xMin, xMax);
   RooRealVar mean("mean", "mean", meanStartVal, xMin, xMax);
   RooRealVar sigma("sigma", "sigma", sigmaStartVal, 0.01, xMax);
   RooGaussian model("model", "model", x, mean, sigma);

   RooArgSet obsSet{x};

   auto dataSetView = RooDataSet::fromArrays("dataSetView", "dataSetView", obsSet, int(nEntries), std::vector<double *>{xValues.data()});

   RooDataSet dataSetCopy{"dataSetCopy", "dataSetCopy", obsSet};
   for (std::size_t i = 0; i < nEntries; ++i) {
      x.setVal(xValues[i]);
      dataSetCopy.add(obsSet);
   }

   std::unique_ptr<RooFitResult> result{
      model.fitTo(*dataSetView, Save(), PrintLevel(-1), Verbose(false), BatchMode(false))};
   mean.setVal(meanStartVal);
   sigma.setVal(sigmaStartVal);
   std::unique_ptr<RooFitResult> resultRef{
      model.fitTo(dataSetCopy, Save(), PrintLevel(-1), Verbose(false), BatchMode(false))};

   EXPECT_TRUE(result->isIdentical(*resultRef));
}
