#include <TFile.h>
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooGlobalFunc.h>

#include <memory>

#include "gtest/gtest.h"

class LoadOldWorkspace : public testing::TestWithParam<std::string> {};

TEST_P(LoadOldWorkspace, DifferentVersions)
{
   TFile file(GetParam().c_str());
   ASSERT_TRUE(file.IsOpen());

   RooWorkspace *w = nullptr;
   file.GetObject("w", w);
   ASSERT_NE(w, nullptr);

   RooAddPdf *model = dynamic_cast<RooAddPdf *>(w->pdf("model"));
   ASSERT_NE(model, nullptr);
   EXPECT_STREQ(model->GetName(), "model");
   EXPECT_STREQ(model->GetTitle(), "g1+g2+a");

   RooDataSet *data = dynamic_cast<RooDataSet *>(w->data("modelData"));
   ASSERT_NE(data, nullptr);

   std::unique_ptr<RooArgSet> observables(model->getObservables(*data));
   std::unique_ptr<RooArgSet> parameters(model->getParameters(*data));

   *observables = *data->get(0);
   EXPECT_NEAR(model->getVal(*observables), 0.393976, 1.E-6);

   *observables = *data->get(1);
   EXPECT_NEAR(model->getVal(*observables), 0.344877, 1.E-6);
}

INSTANTIATE_TEST_SUITE_P(ROOT6, LoadOldWorkspace,
                         testing::Values("rf502_workspace_v6.14.root", "rf502_workspace_v6.04.root"));

INSTANTIATE_TEST_SUITE_P(DISABLED_ROOT5, LoadOldWorkspace, testing::Values("rf502_workspace_v5.34.root"));

namespace {

// Helper to fill a dataset with a few entries.
std::unique_ptr<RooDataSet> makeData(const char *name, RooRealVar &x, std::vector<double> const &vals)
{
   auto data = std::make_unique<RooDataSet>(name, name, x);
   for (double v : vals) {
      x.setVal(v);
      data->add(x);
   }
   return data;
}

} // namespace

// Round-trip test for the current RooWorkspace format, where the owned datasets
// are stored in std::vectors of unique pointers instead of RooLinkedLists (see
// the schema-evolution rules in the LinkDef). This checks that datasets survive
// a write/read cycle, that lookup and insertion order are preserved, and that
// embedded datasets are kept separate.
TEST(RoundTripWorkspace, DataListModernization)
{
   RooRealVar x("x", "x", 0, 10);

   auto d1 = makeData("d1", x, {1., 2.});
   auto d2 = makeData("d2", x, {3., 4., 5.});
   auto emb = makeData("emb", x, {6.});

   RooWorkspace wIn("w");
   ASSERT_FALSE(wIn.import(*d1));
   ASSERT_FALSE(wIn.import(*d2));
   ASSERT_FALSE(wIn.import(*emb, RooFit::Embedded()));

   // Importing a dataset with an already-existing name must be rejected (returns
   // true) and must not change the workspace contents.
   ASSERT_TRUE(wIn.import(*d1));
   EXPECT_EQ(wIn.allData().size(), 2u);

   // A lookup for a non-existing dataset returns a null pointer (and must not crash).
   EXPECT_EQ(wIn.data("does_not_exist"), nullptr);

   // Write and read back.
   TString fname{"roundtrip_workspace.root"};
   {
      TFile fOut(fname, "RECREATE");
      wIn.Write("w");
   }

   TFile fIn(fname);
   ASSERT_TRUE(fIn.IsOpen());
   RooWorkspace *w = nullptr;
   fIn.GetObject("w", w);
   ASSERT_NE(w, nullptr);

   // Regular datasets are found, embedded dataset is not in the regular list.
   ASSERT_NE(w->data("d1"), nullptr);
   ASSERT_NE(w->data("d2"), nullptr);
   EXPECT_EQ(w->data("emb"), nullptr);

   EXPECT_EQ(w->data("d1")->numEntries(), 2);
   EXPECT_EQ(w->data("d2")->numEntries(), 3);

   // Embedded dataset is retrievable through the dedicated accessor.
   ASSERT_NE(w->embeddedData("emb"), nullptr);
   EXPECT_EQ(w->embeddedData("emb")->numEntries(), 1);
   EXPECT_EQ(w->data("emb"), nullptr);

   // Insertion order of the regular datasets is preserved.
   std::list<RooAbsData *> all = w->allData();
   ASSERT_EQ(all.size(), 2u);
   EXPECT_STREQ(all.front()->GetName(), "d1");
   EXPECT_STREQ(all.back()->GetName(), "d2");
}
