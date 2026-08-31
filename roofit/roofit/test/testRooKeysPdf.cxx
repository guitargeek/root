// Tests for the RooKeysPdf and friends
// Authors: Jonas Rembser, CERN  07/2022

#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooKeysPdf.h>
#include <RooNDKeysPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <TMath.h>
#include <TRandom3.h>

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

/// A list of (value, weight) pairs.
using Entries = std::vector<std::pair<double, double>>;

std::unique_ptr<RooDataSet> makeWeightedData(const char *name, RooRealVar &x, RooRealVar &w, Entries const &entries)
{
   auto data = std::make_unique<RooDataSet>(name, name, RooArgSet{x, w}, RooFit::WeightVar(w));
   for (auto const &entry : entries) {
      x.setVal(entry.first);
      data->add(RooArgSet{x, w}, entry.second);
   }
   return data;
}

/// Independent implementation of the kernel estimation that RooKeysPdf does,
/// for the case without mirroring. It follows the documented algorithm:
///
///   * the normal-reference width is `h = (4/3)^(1/5) * W^(-1/5) * rho`, where
///     `W` is the sum of event weights,
///   * the width of the kernel of event `j` is scaled with the inverse square
///     root of the local event density at `x_j`, which is estimated with the
///     *absolute* value of the event weights,
///   * the kernels are summed with the signed event weights and the sum is
///     normalized by `W`.
///
/// The result is sampled on the same grid of `nPoints` values that RooKeysPdf
/// uses internally, and linearly interpolated in between.
class ReferenceKeysPdf {
public:
   ReferenceKeysPdf(Entries const &entries, double lo, double hi, double rho = 1.0) : _lo{lo}, _hi{hi}
   {
      const double sqrt2pi = std::sqrt(2. * TMath::Pi());

      double sumW = 0.;
      double sumWX = 0.;
      double sumWXX = 0.;
      for (auto const &entry : entries) {
         sumW += entry.second;
         sumWX += entry.second * entry.first;
         sumWXX += entry.second * entry.first * entry.first;
      }
      const double mean = sumWX / sumW;
      const double sigma = std::sqrt(sumWXX / sumW - mean * mean);
      const double h = std::pow(4. / 3., 0.2) * std::pow(sumW, -0.2) * rho;
      const double hmin = h * sigma * std::sqrt(2.) / 10.;
      const double norm = h * std::sqrt(sigma * sumW) / (2. * std::sqrt(3.));

      // adaptive kernel widths from the local density of events
      std::vector<double> widths;
      for (auto const &entry : entries) {
         double density = 0.;
         for (auto const &other : entries) {
            const double r = (entry.first - other.first) / (h * sigma);
            density += std::abs(other.second) * std::exp(-0.5 * r * r);
         }
         density /= h * sigma * sqrt2pi;
         widths.push_back(std::max(hmin, norm / std::sqrt(density)));
      }

      // superposition of the kernels
      _table.resize(nPoints + 1, 0.);
      for (std::size_t i = 0; i < entries.size(); ++i) {
         for (int k = 0; k <= nPoints; ++k) {
            const double chi = (_lo + k * binWidth() - entries[i].first) / widths[i] / std::sqrt(2.);
            _table[k] += entries[i].second / widths[i] * std::exp(-chi * chi);
         }
      }
      for (auto &val : _table) {
         val /= sqrt2pi * sumW;
      }
   }

   /// Unnormalized pdf value, clipped at zero like RooKeysPdf::evaluate() does.
   double operator()(double x) const
   {
      int i = std::max(0, std::min(nPoints - 1, static_cast<int>(std::floor((x - _lo) / binWidth()))));
      const double dx = (x - (_lo + i * binWidth())) / binWidth();
      return std::max(0., _table[i] + dx * (_table[i + 1] - _table[i]));
   }

   double max() const { return *std::max_element(_table.begin(), _table.end()); }

private:
   double binWidth() const { return (_hi - _lo) / (nPoints - 1); }

   static constexpr int nPoints = 1000;
   double _lo;
   double _hi;
   std::vector<double> _table;
};

/// Compare a RooKeysPdf to the reference implementation above. The tolerance is
/// given relative to the maximum of the reference density; it has to
/// accommodate that the two implementations sample the kernel sum on slightly
/// different grids before interpolating.
void compareToReference(RooKeysPdf &pdf, RooRealVar &x, Entries const &entries, double relTol, double rho = 1.0)
{
   ReferenceKeysPdf ref{entries, x.getMin(), x.getMax(), rho};
   const double tol = relTol * ref.max();
   for (int i = 0; i <= 200; ++i) {
      const double xVal = x.getMin() + (x.getMax() - x.getMin()) * i / 200.;
      x.setVal(xVal);
      EXPECT_NEAR(pdf.getVal(), ref(xVal), tol) << "at x = " << xVal;
   }
}

} // namespace

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

// Check that the kernel estimation of a dataset with (positive) event weights
// follows the documented algorithm, in particular that the local event density
// that sets the adaptive kernel widths is the weighted one.
TEST(RooKeysPdf, WeightedDatasetReference)
{
   RooRealVar x{"x", "x", 0., 20.};
   RooRealVar w{"w", "w", -1e9, 1e9};

   TRandom3 rng{1234};

   // Mildly varying weights, and weights that differ by two orders of
   // magnitude from event to event.
   Entries mild;
   Entries extreme;
   for (int i = 0; i < 200; ++i) {
      mild.emplace_back(rng.Gaus(10., 1.5), 0.5 + rng.Uniform());
      extreme.emplace_back(rng.Gaus(10., 1.5), i % 2 ? 0.1 : 10.0);
   }

   std::unique_ptr<RooDataSet> dataMild{makeWeightedData("dataMild", x, w, mild)};
   std::unique_ptr<RooDataSet> dataExtreme{makeWeightedData("dataExtreme", x, w, extreme)};

   RooKeysPdf pdfMild{"pdfMild", "", x, *dataMild, RooKeysPdf::NoMirror};
   RooKeysPdf pdfExtreme{"pdfExtreme", "", x, *dataExtreme, RooKeysPdf::NoMirror};

   compareToReference(pdfMild, x, mild, 0.02);
   compareToReference(pdfExtreme, x, extreme, 0.03);
}

// Datasets with negative event weights are supported, as long as the estimation
// is well defined. This covers GitHub issue #12639.
TEST(RooKeysPdf, NegativeWeights)
{
   RooRealVar x{"x", "x", 0., 20.};
   RooRealVar w{"w", "w", -1e9, 1e9};

   TRandom3 rng{4321};

   // A dataset that mimics background subtraction with per-event weights, as
   // done with sWeights: a signal sample with positive weights, plus a
   // background sample that is statistically cancelled by an independent
   // background sample with negative weights. In expectation, the weighted
   // density is the signal density.
   const int nSig = 2000;
   const int nBkg = 1000;
   Entries subtracted;
   Entries unsubtracted;
   for (int i = 0; i < nSig; ++i) {
      const double val = rng.Gaus(10., 1.5);
      subtracted.emplace_back(val, 1.0);
      unsubtracted.emplace_back(val, 1.0);
   }
   for (int i = 0; i < nBkg; ++i) {
      const double val = rng.Uniform(0., 20.);
      subtracted.emplace_back(val, 1.0);
      unsubtracted.emplace_back(val, 1.0);
   }
   for (int i = 0; i < nBkg; ++i) {
      subtracted.emplace_back(rng.Uniform(0., 20.), -1.0);
   }

   std::unique_ptr<RooDataSet> data{makeWeightedData("dataSub", x, w, subtracted)};
   std::unique_ptr<RooDataSet> dataUnsub{makeWeightedData("dataUnsub", x, w, unsubtracted)};

   RooKeysPdf pdf{"pdfSub", "", x, *data, RooKeysPdf::NoMirror};
   RooKeysPdf pdfUnsub{"pdfUnsub", "", x, *dataUnsub, RooKeysPdf::NoMirror};

   RooArgSet normSet{x};

   // The estimation must be finite and non-negative everywhere, and it has to
   // follow the same algorithm as for positive weights.
   for (int i = 0; i <= 200; ++i) {
      x.setVal(x.getMin() + (x.getMax() - x.getMin()) * i / 200.);
      const double val = pdf.getVal(normSet);
      EXPECT_TRUE(std::isfinite(val)) << "at x = " << x.getVal();
      EXPECT_GE(val, 0.0) << "at x = " << x.getVal();
   }
   compareToReference(pdf, x, subtracted, 0.03);

   // What the user gets out of the pdf has to integrate to unity. This only
   // works if the clipping of the negative parts of the estimation is applied
   // consistently in RooKeysPdf::evaluate() and in the analytical integral that
   // normalizes it.
   {
      const int nSteps = 20001;
      const double step = (x.getMax() - x.getMin()) / nSteps;
      double norm = 0.;
      for (int i = 0; i <= nSteps; ++i) {
         x.setVal(x.getMin() + i * step);
         norm += pdf.getVal(normSet) * step * (i == 0 || i == nSteps ? 0.5 : 1.0);
      }
      EXPECT_NEAR(norm, 1.0, 1e-3);
   }

   // The subtracted estimation has to be much closer to the true signal shape
   // than the unsubtracted one.
   double distSub = 0.;
   double distUnsub = 0.;
   const double step = (x.getMax() - x.getMin()) / 1000.;
   for (int i = 0; i <= 1000; ++i) {
      const double xVal = x.getMin() + i * step;
      x.setVal(xVal);
      const double r = (xVal - 10.) / 1.5;
      const double truth = std::exp(-0.5 * r * r) / (1.5 * std::sqrt(2. * TMath::Pi()));
      distSub += std::abs(pdf.getVal(normSet) - truth) * step;
      distUnsub += std::abs(pdfUnsub.getVal(normSet) - truth) * step;
   }
   EXPECT_LT(distSub, 0.15);
   EXPECT_GT(distUnsub, 3. * distSub);
}

// The event weights of a RooKeysPdf input dataset are event multiplicities: an
// event with weight two has to be exactly equivalent to two events of weight one
// at the same value. This is the convention that RooKeysPdf follows since the
// support for weighted datasets was introduced, and the RooKeysPdf.WeightedDataset
// test above relies on it. Before the fix for GitHub issue #12639 it only held
// when all events in a neighbourhood happened to carry the same weight, because
// the local density that sets the kernel widths was the event's own weight times
// the *unweighted* density instead of the weighted one.
TEST(RooKeysPdf, WeightsAreMultiplicities)
{
   RooRealVar x{"x", "x", 0., 20.};
   RooRealVar w{"w", "w", -1e9, 1e9};

   TRandom3 rng{9876};

   Entries compact;
   Entries expanded;
   for (int i = 0; i < 80; ++i) {
      double val = 0.;
      do {
         val = rng.Gaus(10., 2.0);
      } while (val <= x.getMin() || val >= x.getMax());
      const int multiplicity = 1 + i % 3;
      compact.emplace_back(val, double(multiplicity));
      for (int k = 0; k < multiplicity; ++k) {
         expanded.emplace_back(val, 1.0);
      }
   }

   std::unique_ptr<RooDataSet> dataCompact{makeWeightedData("dataCompact", x, w, compact)};
   std::unique_ptr<RooDataSet> dataExpanded{makeWeightedData("dataExpanded", x, w, expanded)};

   RooKeysPdf pdfCompact{"pdfCompact", "", x, *dataCompact, RooKeysPdf::NoMirror};
   RooKeysPdf pdfExpanded{"pdfExpanded", "", x, *dataExpanded, RooKeysPdf::NoMirror};
   RooNDKeysPdf ndCompact{"ndCompact", "", x, *dataCompact, "a"};
   RooNDKeysPdf ndExpanded{"ndExpanded", "", x, *dataExpanded, "a"};

   RooArgSet normSet{x};
   for (int i = 0; i <= 200; ++i) {
      x.setVal(x.getMin() + (x.getMax() - x.getMin()) * i / 200.);
      EXPECT_NEAR(pdfCompact.getVal(normSet), pdfExpanded.getVal(normSet), 1e-8) << "at x = " << x.getVal();
      EXPECT_NEAR(ndCompact.getVal(normSet), ndExpanded.getVal(normSet), 1e-8) << "at x = " << x.getVal();
   }
}

// Negative event weights can make the kernel estimation ill-defined. Such cases
// have to be rejected with an explanatory exception instead of silently
// producing NaN pdf values. This covers GitHub issue #12639.
TEST(RooKeysPdf, ThrowOnIllDefinedInput)
{
   RooRealVar x{"x", "x", 0., 20.};
   RooRealVar w{"w", "w", -1e9, 1e9};

   TRandom3 rng{1};

   Entries zeroSumW;
   Entries negativeSumW;
   for (int i = 0; i < 50; ++i) {
      zeroSumW.emplace_back(rng.Gaus(10., 2.), 1.0);
      zeroSumW.emplace_back(rng.Gaus(10., 2.), -1.0);
      negativeSumW.emplace_back(rng.Gaus(10., 2.), 1.0);
      negativeSumW.emplace_back(rng.Gaus(10., 2.), -2.0);
   }
   // Positive total weight, but the negative weights make the weighted variance
   // negative.
   Entries negativeVariance{{10.0, 10.0}, {10.1, 10.0}, {2.0, -1.0}, {18.0, -1.0}};
   // All events at the same value: no spread to estimate a kernel width from.
   Entries noSpread{{10.0, 1.0}, {10.0, 2.0}, {10.0, 3.0}};

   std::unique_ptr<RooDataSet> dataZero{makeWeightedData("dataZero", x, w, zeroSumW)};
   std::unique_ptr<RooDataSet> dataNeg{makeWeightedData("dataNeg", x, w, negativeSumW)};
   std::unique_ptr<RooDataSet> dataVar{makeWeightedData("dataVar", x, w, negativeVariance)};
   std::unique_ptr<RooDataSet> dataFlat{makeWeightedData("dataFlat", x, w, noSpread)};
   RooDataSet dataEmpty{"dataEmpty", "dataEmpty", x};

   auto errorMessage = [&](RooDataSet &data) {
      try {
         RooKeysPdf pdf{"pdf", "", x, data, RooKeysPdf::NoMirror};
      } catch (std::runtime_error const &e) {
         return std::string{e.what()};
      }
      return std::string{};
   };

   EXPECT_NE(errorMessage(*dataZero).find("total weight"), std::string::npos);
   EXPECT_NE(errorMessage(*dataNeg).find("total weight"), std::string::npos);
   EXPECT_NE(errorMessage(dataEmpty).find("total weight"), std::string::npos);
   EXPECT_NE(errorMessage(*dataVar).find("variance"), std::string::npos);
   EXPECT_NE(errorMessage(*dataFlat).find("variance"), std::string::npos);

   // The same has to hold for the multi-dimensional version.
   EXPECT_THROW(RooNDKeysPdf("ndZero", "", x, *dataZero, "a"), std::runtime_error);
   EXPECT_THROW(RooNDKeysPdf("ndNeg", "", x, *dataNeg, "a"), std::runtime_error);
   EXPECT_THROW(RooNDKeysPdf("ndVar", "", x, *dataVar, "a"), std::runtime_error);
   EXPECT_THROW(RooNDKeysPdf("ndFlat", "", x, *dataFlat, "a"), std::runtime_error);
}

// RooNDKeysPdf has to cope with negative event weights as well.
TEST(RooNDKeysPdf, NegativeWeights)
{
   RooRealVar x{"x", "x", 0., 20.};
   RooRealVar w{"w", "w", -1e9, 1e9};

   TRandom3 rng{4321};

   Entries entries;
   for (int i = 0; i < 2000; ++i) {
      entries.emplace_back(rng.Gaus(10., 1.5), 1.0);
   }
   for (int i = 0; i < 1000; ++i) {
      entries.emplace_back(rng.Uniform(0., 20.), 1.0);
      entries.emplace_back(rng.Uniform(0., 20.), -1.0);
   }

   std::unique_ptr<RooDataSet> data{makeWeightedData("dataNDSub", x, w, entries)};
   RooNDKeysPdf pdf{"pdfND", "", x, *data, "a"};

   RooArgSet normSet{x};
   double dist = 0.;
   const double step = (x.getMax() - x.getMin()) / 1000.;
   for (int i = 0; i <= 1000; ++i) {
      const double xVal = x.getMin() + i * step;
      x.setVal(xVal);
      const double val = pdf.getVal(normSet);
      EXPECT_TRUE(std::isfinite(val)) << "at x = " << xVal;
      EXPECT_GE(val, 0.0) << "at x = " << xVal;
      const double r = (xVal - 10.) / 1.5;
      dist += std::abs(val - std::exp(-0.5 * r * r) / (1.5 * std::sqrt(2. * TMath::Pi()))) * step;
   }
   EXPECT_LT(dist, 0.15);
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
