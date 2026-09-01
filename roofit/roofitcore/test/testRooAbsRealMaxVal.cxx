// Tests for the RooAbsReal::getMaxVal() / RooAbsReal::maxVal() contract,
// see https://github.com/root-project/root/issues/12317.
//
// The contract (documented in RooAbsReal::getMaxVal()) is that
// maxVal(getMaxVal(vars)) is an upper bound on the *unnormalized* function
// value while the observables in `vars` scan their full range. The tests below
// check that property by densely scanning the observables.
//
// Author: Jonas Rembser, CERN

#include <RooBinSamplingPdf.h>
#include <RooBinning.h>
#include <RooCBShape.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFFTConvPdf.h>
#include <RooGaussian.h>
#include <RooHistFunc.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooNumGenConfig.h>
#include <RooRandom.h>
#include <RooRealProxy.h>
#include <RooRealVar.h>
#include <RooStudentT.h>
#include <RooWrapperPdf.h>

#ifdef ROOFITMORE
#include <RooLegendre.h>
#include <RooSpHarmonic.h>
#endif

#include <TMath.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

/// Maximum of the unnormalized function value, scanning `vars` on a dense grid.
/// For functions that can be negative, the maximum of the absolute value is
/// what maxVal() bounds.
double scanMax(RooAbsReal const &func, std::vector<RooRealVar *> const &vars, int nPoints, bool useAbs)
{
   RooArgSet emptyNormSet{};
   double max = -std::numeric_limits<double>::infinity();

   std::vector<int> idx(vars.size(), 0);
   while (true) {
      for (std::size_t i = 0; i < vars.size(); ++i) {
         double lo = vars[i]->getMin();
         double hi = vars[i]->getMax();
         vars[i]->setVal(lo + (hi - lo) * idx[i] / double(nPoints));
      }
      const double val = func.getVal(emptyNormSet);
      max = std::max(max, useAbs ? std::abs(val) : val);

      // Increment the odometer.
      std::size_t dim = 0;
      for (; dim < vars.size(); ++dim) {
         if (++idx[dim] <= nPoints)
            break;
         idx[dim] = 0;
      }
      if (dim == vars.size())
         break;
   }
   return max;
}

/// Assert that the advertised maximum is an upper bound on the function value,
/// and that it is not looser than `maxOverestimation` times the true maximum.
/// Note that the scan can only get close to the true maximum, so even an exact
/// bound needs a small tolerance in `maxOverestimation`.
void checkMaxVal(RooAbsReal const &func, std::vector<RooRealVar *> const &vars, double maxOverestimation,
                 int nPoints = 10001, bool useAbs = false)
{
   RooArgSet varSet{};
   for (RooRealVar *var : vars) {
      varSet.add(*var);
   }

   const int code = func.getMaxVal(varSet);
   ASSERT_NE(code, 0) << func.ClassName() << "(" << func.GetName() << ") does not advertise a maximum";
   const double maxVal = func.maxVal(code);

   const double trueMax = scanMax(func, vars, nPoints, useAbs);

   // Underestimating is a correctness bug: it silently biases accept-reject sampling.
   EXPECT_GE(maxVal, trueMax) << func.ClassName() << "(" << func.GetName() << ") underestimates its maximum";
   // Overestimating only costs generation efficiency, but should stay bounded.
   EXPECT_LE(maxVal, maxOverestimation * trueMax)
      << func.ClassName() << "(" << func.GetName() << ") overestimates its maximum too much";
}

/// The maximum must not be advertised for an empty set of scanned observables,
/// and it must be advertised (and unchanged) if the set contains additional
/// observables that the function does not depend on.
void checkNormSetHandling(RooAbsReal const &func, RooArgSet const &obs, RooAbsArg &disconnected)
{
   EXPECT_EQ(func.getMaxVal(RooArgSet{}), 0)
      << func.ClassName() << "(" << func.GetName() << ") advertises a maximum for an empty set of observables";

   const int code = func.getMaxVal(obs);
   ASSERT_NE(code, 0);

   RooArgSet obsPlusExtra{obs};
   obsPlusExtra.add(disconnected);
   const int codeExtra = func.getMaxVal(obsPlusExtra);
   ASSERT_NE(codeExtra, 0) << func.ClassName() << "(" << func.GetName()
                           << ") does not handle extra disconnected observables";
   EXPECT_DOUBLE_EQ(func.maxVal(codeExtra), func.maxVal(code));
}

/// A pdf that deliberately advertises a maximum that is too small, used to
/// check that the accept-reject generator complains about it.
class LyingPdf : public RooAbsPdf {
public:
   LyingPdf(const char *name, RooAbsReal &x) : RooAbsPdf{name, name}, _x{"x", "x", this, x} {}
   LyingPdf(const LyingPdf &other, const char *name = nullptr) : RooAbsPdf{other, name}, _x{"x", this, other._x} {}
   TObject *clone(const char *newname = nullptr) const override { return new LyingPdf{*this, newname}; }

   Int_t getMaxVal(const RooArgSet &vars) const override { return vars.contains(*_x.absArg()) ? 1 : 0; }
   /// The true maximum is one, reached at x == 0.5.
   double maxVal(Int_t) const override { return 0.1; }

protected:
   double evaluate() const override { return std::exp(-0.5 * std::pow((_x - 0.5) / 0.2, 2)); }

private:
   RooRealProxy _x;
};

/// Make the accept-reject sampler the one-dimensional generator. The category
/// states of the generator configuration are only defined once the
/// RooNumGenFactory has been instantiated, which happens on the first
/// generation, hence the warm-up.
bool selectAcceptRejectSampler()
{
   // A RooCBShape has no internal generator, so this goes through the numeric
   // sampler factory.
   RooRealVar x{"x", "x", 0.5, 0., 1.};
   RooRealVar mean{"mean", "mean", 0.5};
   RooRealVar sigma{"sigma", "sigma", 0.2};
   RooRealVar alpha{"alpha", "alpha", 1.0};
   RooRealVar n{"n", "n", 1.0};
   RooCBShape cbShape{"cbShape", "cbShape", x, mean, sigma, alpha, n};
   std::unique_ptr<RooDataSet> warmUp{cbShape.generate(x, 1)};

   RooNumGenConfig::defaultConfig().method1D(false, false).setLabel("RooAcceptReject");
   return std::string("RooAcceptReject") == RooNumGenConfig::defaultConfig().method1D(false, false).getCurrentLabel();
}

} // namespace

/// The crystal ball shapes are unimodal with a maximum of exactly one. Before
/// the fix for issue #12317, they returned one over the normalization
/// integral, which underestimates the maximum whenever that integral is larger
/// than one.
TEST(RooAbsRealMaxVal, CrystalBall)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   RooRealVar y{"y", "y", 2.5, 0, 5};
   RooRealVar x0{"x0", "x0", 2.5};
   RooRealVar sigma{"sigma", "sigma", 0.02};
   RooRealVar alpha{"alpha", "alpha", 1.0};
   RooRealVar n{"n", "n", 1.0};

   RooCrystalBall crystalBall{"crystalBall", "crystalBall", x, x0, sigma, alpha, n};
   RooCBShape cbShape{"cbShape", "cbShape", x, x0, sigma, alpha, n};
   RooCrystalBall doubleSided{"doubleSided", "doubleSided", x, x0, sigma, alpha, n, /*doubleSided=*/true};

   // A negative alpha puts the tail on the other side, both for RooCBShape and
   // for the single-tail RooCrystalBall, and values of n at and below one hit
   // the special cases of the power-law tail.
   for (double sigmaVal : {0.02, 0.5, 2.0}) {
      for (double alphaVal : {-5.0, -1.0, -0.3, 0.001, 0.3, 1.0, 5.0}) {
         for (double nVal : {0.0, 1e-6, 0.5, 1.0, 10.0}) {
            sigma.setVal(sigmaVal);
            alpha.setVal(alphaVal);
            n.setVal(nVal);
            checkMaxVal(crystalBall, {&x}, 1.001);
            checkMaxVal(cbShape, {&x}, 1.001);
            if (alphaVal > 0.0) {
               // The two-sided shape requires a positive alpha.
               checkMaxVal(doubleSided, {&x}, 1.001);
            }
         }
      }
   }

   // A crystal ball with different widths and different tails on the two sides.
   RooRealVar sigmaL{"sigmaL", "sigmaL", 0.3};
   RooRealVar sigmaR{"sigmaR", "sigmaR", 1.5};
   RooRealVar alphaL{"alphaL", "alphaL", 0.5};
   RooRealVar nL{"nL", "nL", 0.2};
   RooRealVar alphaR{"alphaR", "alphaR", 2.0};
   RooRealVar nR{"nR", "nR", 5.0};
   RooCrystalBall asymmetric{"asymmetric", "asymmetric", x, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR};
   for (double sigmaLVal : {0.05, 0.3, 2.0}) {
      for (double alphaLVal : {0.05, 0.5, 3.0}) {
         for (double nLVal : {0.0, 0.2, 1.0, 7.0}) {
            sigmaL.setVal(sigmaLVal);
            alphaL.setVal(alphaLVal);
            nL.setVal(nLVal);
            checkMaxVal(asymmetric, {&x}, 1.001);
         }
      }
   }

   sigma.setVal(0.02);
   alpha.setVal(1.0);
   n.setVal(1.0);
   checkNormSetHandling(crystalBall, RooArgSet{x}, y);
   checkNormSetHandling(cbShape, RooArgSet{x}, y);
   checkNormSetHandling(doubleSided, RooArgSet{x}, y);
}

TEST(RooAbsRealMaxVal, StudentT)
{
   RooRealVar x{"x", "x", 0.0, -10, 10};
   RooRealVar y{"y", "y", 0.0, -10, 10};
   RooRealVar mean{"mean", "mean", 0.0};
   RooRealVar sigma{"sigma", "sigma", 1.5};
   RooRealVar ndf{"ndf", "ndf", 3.0};
   RooStudentT studentT{"studentT", "studentT", x, mean, sigma, ndf};

   for (double ndfVal : {1.0, 3.0, 10.0}) {
      ndf.setVal(ndfVal);
      checkMaxVal(studentT, {&x}, 1.001);
   }
   checkNormSetHandling(studentT, RooArgSet{x}, y);
}

/// A RooHistPdf divides the bin content by the bin volume, a RooHistFunc does
/// not. Before the fix for issue #12317, both returned the largest bin content,
/// so the RooHistPdf underestimated its maximum by the bin volume.
TEST(RooAbsRealMaxVal, HistPdfAndHistFunc)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   RooRealVar y{"y", "y", 2.5, 0, 5};
   x.setBins(25);

   RooDataHist dataHist{"dataHist", "dataHist", x};
   for (int i = 0; i < x.numBins(); ++i) {
      x.setBin(i);
      dataHist.set(x, 10. + 90. * std::exp(-0.5 * std::pow((x.getVal() - 2.5) / 0.5, 2)), -1.);
   }
   x.setVal(2.5);

   // The 1.05 safety margin that the implementation adds explains the tolerance.
   for (int intOrder : {0, 1}) {
      RooHistFunc histFunc{"histFunc", "histFunc", x, dataHist, intOrder};
      RooHistPdf histPdf{"histPdf", "histPdf", x, dataHist, intOrder};
      RooWrapperPdf wrapperPdf{"wrapperPdf", "wrapperPdf", histFunc};

      checkMaxVal(histFunc, {&x}, 1.06);
      checkMaxVal(histPdf, {&x}, 1.06);
      checkMaxVal(wrapperPdf, {&x}, 1.06);

      checkNormSetHandling(histFunc, RooArgSet{x}, y);
      checkNormSetHandling(histPdf, RooArgSet{x}, y);
      checkNormSetHandling(wrapperPdf, RooArgSet{x}, y);
   }

   // Interpolation orders higher than one are not bounded by the bin contents,
   // so no maximum must be advertised.
   {
      RooHistFunc histFunc{"histFunc2", "histFunc2", x, dataHist, 2};
      RooHistPdf histPdf{"histPdf2", "histPdf2", x, dataHist, 2};
      EXPECT_EQ(histFunc.getMaxVal(RooArgSet{x}), 0);
      EXPECT_EQ(histPdf.getMaxVal(RooArgSet{x}), 0);
   }

   // The c.d.f. boundary conditions extend the histogram with a bin of content
   // one, so with interpolation the bin contents are not a bound any more.
   {
      RooHistFunc histFunc{"histFunc3", "histFunc3", x, dataHist, 1};
      RooHistPdf histPdf{"histPdf3", "histPdf3", x, dataHist, 1};
      histFunc.setCdfBoundaries(true);
      histPdf.setCdfBoundaries(true);
      EXPECT_EQ(histFunc.getMaxVal(RooArgSet{x}), 0);
      EXPECT_EQ(histPdf.getMaxVal(RooArgSet{x}), 0);

      // Without interpolation the boundary conditions play no role.
      RooHistFunc histFunc0{"histFunc4", "histFunc4", x, dataHist, 0};
      histFunc0.setCdfBoundaries(true);
      checkMaxVal(histFunc0, {&x}, 1.06);
   }
}

/// The bound must stay a bound if the bin contents are negative, where scaling
/// it up by the safety margin would push it below the true maximum. A
/// RooHistPdf clips its value at zero, a RooHistFunc does not.
TEST(RooAbsRealMaxVal, HistPdfAndHistFuncNegativeBins)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   x.setBins(5);

   RooDataHist dataHist{"dataHist", "dataHist", x};
   for (int i = 0; i < x.numBins(); ++i) {
      x.setBin(i);
      dataHist.set(x, -10. - i, -1.);
   }
   x.setVal(2.5);

   RooHistFunc histFunc{"histFunc", "histFunc", x, dataHist, 0};
   RooHistPdf histPdf{"histPdf", "histPdf", x, dataHist, 0};

   const double funcMax = histFunc.maxVal(histFunc.getMaxVal(RooArgSet{x}));
   const double pdfMax = histPdf.maxVal(histPdf.getMaxVal(RooArgSet{x}));

   RooArgSet emptyNormSet{};
   for (int i = 0; i < x.numBins(); ++i) {
      x.setBin(i);
      EXPECT_LE(histFunc.getVal(emptyNormSet), funcMax) << "RooHistFunc underestimates its maximum in bin " << i;
      EXPECT_LE(histPdf.getVal(emptyNormSet), pdfMax) << "RooHistPdf underestimates its maximum in bin " << i;
   }
}

/// With non-uniform binning, the largest bin content is not the largest
/// density, which is what a RooHistPdf evaluates to.
TEST(RooAbsRealMaxVal, HistPdfNonUniformBinning)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   const double edges[6] = {0., 0.1, 0.5, 2.0, 3.0, 5.0};
   x.setBinning(RooBinning(5, edges));

   RooDataHist dataHist{"dataHist", "dataHist", x};
   for (int i = 0; i < x.numBins(); ++i) {
      x.setBin(i);
      dataHist.set(x, 10. * (i + 1), -1.);
   }
   x.setVal(2.5);

   RooHistFunc histFunc{"histFunc", "histFunc", x, dataHist, 0};
   RooHistPdf histPdf{"histPdf", "histPdf", x, dataHist, 0};

   checkMaxVal(histFunc, {&x}, 1.06);
   checkMaxVal(histPdf, {&x}, 1.06);
}

TEST(RooAbsRealMaxVal, KeysPdf)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   RooRealVar y{"y", "y", 2.5, 0, 5};

   RooRandom::randomGenerator()->SetSeed(1337);
   RooDataSet data{"data", "data", x};
   for (int i = 0; i < 500; ++i) {
      x.setVal(std::min(4.999, std::max(0.001, RooRandom::randomGenerator()->Gaus(2.5, 0.6))));
      data.add(x);
   }
   x.setVal(2.5);

   RooKeysPdf keysPdf{"keysPdf", "keysPdf", x, data};
   checkMaxVal(keysPdf, {&x}, 1.001);
   checkNormSetHandling(keysPdf, RooArgSet{x}, y);

   // If the observable of the pdf can leave the range that the internal lookup
   // table was built for, evaluate() extrapolates and the table maximum is not
   // a bound anymore, so no maximum may be advertised.
   RooRealVar xWide{"x", "x", 2.5, -20, 25};
   RooKeysPdf keysPdfWide{"keysPdfWide", "keysPdfWide", xWide, x, data};
   EXPECT_EQ(keysPdfWide.getMaxVal(RooArgSet{xWide}), 0);
}

/// The value of a RooBinSamplingPdf is the wrapped pdf averaged over a bin,
/// which is not bounded by the maximum of the wrapped pdf in a way that this
/// class could reliably exploit, so it must not advertise a maximum (see
/// issue #12317).
TEST(RooAbsRealMaxVal, BinSamplingPdfDoesNotAdvertise)
{
   RooRealVar x{"x", "x", 2.5, 0, 5};
   x.setBins(25);
   RooRealVar x0{"x0", "x0", 2.5};
   RooRealVar sigma{"sigma", "sigma", 2.0};
   RooRealVar alpha{"alpha", "alpha", 1.0};
   RooRealVar n{"n", "n", 1.0};
   RooCrystalBall crystalBall{"crystalBall", "crystalBall", x, x0, sigma, alpha, n};
   RooBinSamplingPdf binSamplingPdf{"binSamplingPdf", "binSamplingPdf", x, crystalBall};

   EXPECT_EQ(binSamplingPdf.getMaxVal(RooArgSet{x}), 0);
}

/// A convolution is not bounded by the maximum of the first input pdf, so
/// RooFFTConvPdf must not forward it (see issue #12317).
TEST(RooAbsRealMaxVal, FFTConvPdfDoesNotAdvertise)
{
   RooRealVar x{"x", "x", 0.0, -10, 10};
   RooRealVar m1{"m1", "m1", 0.0};
   RooRealVar s1{"s1", "s1", 1.0};
   RooRealVar a1{"a1", "a1", 1.0};
   RooRealVar n1{"n1", "n1", 2.0};
   // RooCBShape does implement getMaxVal(), so before the fix this was forwarded.
   RooCBShape phys{"phys", "phys", x, m1, s1, a1, n1};
   RooRealVar m2{"m2", "m2", 0.0};
   RooRealVar s2{"s2", "s2", 1.0};
   RooGaussian resolution{"resolution", "resolution", x, m2, s2};
   RooFFTConvPdf fft{"fft", "fft", x, phys, resolution};

   EXPECT_EQ(fft.getMaxVal(RooArgSet{x}), 0);
}

#ifdef ROOFITMORE
TEST(RooAbsRealMaxVal, Legendre)
{
   RooRealVar cosTheta{"cosTheta", "cosTheta", 0.0, -1, 1};
   RooRealVar y{"y", "y", 0.0, -1, 1};

   for (int l = 0; l < 3; ++l) {
      for (int m = 0; m <= l; ++m) {
         RooLegendre legendre{"legendre", "legendre", cosTheta, l, m};
         // The bound for the (l=2, m=1) case is a factor two loose.
         checkMaxVal(legendre, {&cosTheta}, 2.01, 10001, /*useAbs=*/true);
         checkNormSetHandling(legendre, RooArgSet{cosTheta}, y);
      }
   }
}

TEST(RooAbsRealMaxVal, SpHarmonic)
{
   RooRealVar cosTheta{"cosTheta", "cosTheta", 0.0, -1, 1};
   RooRealVar phi{"phi", "phi", 0.0, -TMath::Pi(), TMath::Pi()};

   for (int l = 0; l < 3; ++l) {
      for (int m = 0; m <= l; ++m) {
         RooSpHarmonic spHarmonic{"spHarmonic", "spHarmonic", cosTheta, phi, l, m};
         checkMaxVal(spHarmonic, {&cosTheta, &phi}, 2.01, 501, /*useAbs=*/true);
      }
   }
}
#endif

/// A maximum that is too small silently truncates the generated distribution,
/// because every point above it gets accepted with the same probability. Make
/// sure the accept-reject generator complains loudly instead (see issue #12317).
TEST(RooAbsRealMaxVal, GeneratorComplainsAboutTooSmallMaximum)
{
   const std::string oldMethod = RooNumGenConfig::defaultConfig().method1D(false, false).getCurrentLabel();
   if (!selectAcceptRejectSampler()) {
      GTEST_SKIP() << "could not select the RooAcceptReject sampler";
   }
   struct RestoreMethod {
      ~RestoreMethod() { RooNumGenConfig::defaultConfig().method1D(false, false).setLabel(_label.c_str()); }
      std::string _label;
   } restoreMethod{oldMethod};

   RooRealVar x{"x", "x", 0.5, 0., 1.};
   LyingPdf lyingPdf{"lyingPdf", x};

   RooRandom::randomGenerator()->SetSeed(42);
   testing::internal::CaptureStdout();
   std::unique_ptr<RooDataSet> data{lyingPdf.generate(x, 100)};
   const std::string output = testing::internal::GetCapturedStdout();

   EXPECT_NE(output.find("that was advertised by RooAbsReal::maxVal()"), std::string::npos)
      << "The generator did not complain about the too small maximum. Output was:\n"
      << output;
}
