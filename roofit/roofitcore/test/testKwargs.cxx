// Tests for the overloads of RooFit functions that take kwargs structs
// Authors: Jonas Rembser, CERN 11/2026

#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFit/Kwargs.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHelpers.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <TH1.h>

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

namespace {

// The command arguments of the kwargs struct, by opcode.
template <class Kwargs_t>
std::vector<std::string> opcodes(Kwargs_t const &kwargs)
{
   std::vector<std::string> out;
   for (RooCmdArg const *arg : RooFit::Detail::collectCmdArgs(kwargs, {})) {
      out.emplace_back(arg->opcode());
   }
   return out;
}

} // namespace

// The members of a kwargs struct are command arguments that are empty unless
// given, and that can be initialized like the RooFit function of the same name.
TEST(Kwargs, Members)
{
   RooAbsReal::PlotOnKwargs kwargs;
   EXPECT_TRUE(opcodes(kwargs).empty());

   kwargs.Range = "sideband";
   kwargs.MoveToBack = true;
   EXPECT_EQ(opcodes(kwargs), (std::vector<std::string>{"MoveToBack", "RangeWithName"}));
   EXPECT_STREQ(kwargs.Range.getString(0), "sideband");

   kwargs.Range = {-1., 2.};
   kwargs.MoveToBack = false;
   EXPECT_EQ(opcodes(kwargs), std::vector<std::string>{"Range"});
   EXPECT_EQ(kwargs.Range.getDouble(0), -1.);
   EXPECT_EQ(kwargs.Range.getDouble(1), 2.);

   RooAbsData::PlotOnKwargs dataKwargs;
   dataKwargs.DataError = nullptr;
   EXPECT_EQ(dataKwargs.DataError.getInt(0), RooFit::DataError("None").getInt(0));
}

// Fitting with kwargs, mixed with command arguments, gives the same result as
// with command arguments only.
TEST(Kwargs, FitTo)
{
   RooHelpers::LocalChangeMsgLevel changeMsgLvl(RooFit::WARNING);

   RooRealVar x("x", "x", -10, 10);
   RooRealVar mu("mu", "mu", 1, -10, 10);
   RooRealVar sigma("sigma", "sigma", 2, 0.1, 10);
   RooGaussian gauss("gauss", "gauss", x, mu, sigma);
   x.setRange("sideband", -10, 0);

   std::unique_ptr<RooDataSet> data{gauss.generate(x, 1000)};

   auto fit = [&](auto &&...args) {
      mu.setVal(1.);
      sigma.setVal(2.);
      mu.setError(0.);
      sigma.setError(0.);
      return std::unique_ptr<RooFitResult>{gauss.fitTo(*data, args...)};
   };

   RooFit::FitToKwargs kwargs;
   kwargs.Range = "sideband";
   kwargs.Save = true;
   auto res1 = fit(kwargs, RooFit::PrintLevel(-1));
   auto res2 = fit(RooFit::Range("sideband"), RooFit::Save(), RooFit::PrintLevel(-1));

   ASSERT_NE(res1, nullptr);
   ASSERT_NE(res2, nullptr);
   EXPECT_TRUE(res1->isIdentical(*res2));
}

// Constructors that take kwargs.
TEST(Kwargs, Constructors)
{
   RooRealVar x("x", "x", 0, 10);
   RooRealVar w("w", "w", 0, 10);

   RooDataSet::Kwargs kwargs;
   kwargs.WeightVar = w;
   RooDataSet data{"data", "data", {x, w}, kwargs};
   EXPECT_TRUE(data.isWeighted());

   TH1D hist{"hist", "hist", 10, 0, 10};
   hist.Fill(5., 3.);
   RooDataHist::Kwargs histKwargs;
   histKwargs.Import = hist;
   RooDataHist dataHist{"dataHist", "dataHist", RooArgList{x}, histKwargs};
   EXPECT_DOUBLE_EQ(dataHist.sumEntries(), 3.);
}

// Functions without a RooLinkedList overload take only as many command
// arguments as their overload that takes them one by one.
TEST(Kwargs, TooManyCommandArguments)
{
   RooRealVar x("x", "x", 0, 10);
   RooDataSet data{"data", "data", x};

   RooAbsData::ReduceKwargs kwargs;
   kwargs.Cut = "x > 1";
   kwargs.Name = "reduced";
   RooCmdArg cmd = RooFit::Title("title");
   std::unique_ptr<RooAbsData> reduced{data.reduce(kwargs, cmd)};
   EXPECT_STREQ(reduced->GetName(), "reduced");
   EXPECT_THROW(data.reduce(kwargs, cmd, cmd, cmd, cmd, cmd, cmd, cmd), std::invalid_argument);
}

// The kwargs of the global RooFit functions end up in the sub-arguments of the
// resulting command argument.
TEST(Kwargs, GlobalFunctions)
{
   RooRealVar x("x", "x", 0, 10);
   RooRealVar y("y", "y", 0, 10);

   RooFit::YVarKwargs kwargs;
   kwargs.Binning = 20;
   RooCmdArg yVar = RooFit::YVar(y, kwargs);
   ASSERT_EQ(yVar.subArgs().GetSize(), 1);
   EXPECT_STREQ(yVar.subArgs().At(0)->GetName(), "BinningSpec");

   std::unique_ptr<TH1> hist{x.createHistogram("hist", yVar)};
   EXPECT_EQ(hist->GetNbinsY(), 20);
}
