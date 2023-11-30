#include "gtest/gtest.h"

#include "TH2Poly.h"
#include "TRandom3.h"

TH2Poly *createPoly(Double_t weight = 1)
{
   TH2Poly *h2p = new TH2Poly();
   Int_t i, j;
   Int_t nx = 40;
   Int_t ny = 40;
   Double_t xval1, yval1, xval2, yval2;
   Double_t dx = 0.2, dy = 0.1;
   xval1 = 0.;
   xval2 = dx;
   for (i = 0; i < nx; i++) {
      yval1 = 0.;
      yval2 = dy;
      for (j = 0; j < ny; j++) {
         h2p->AddBin(xval1, yval1, xval2, yval2);
         yval1 = yval2;
         yval2 = yval2 + yval2 * dy;
      }
      xval1 = xval2;
      xval2 = xval2 + xval2 * dx;
   }
   TRandom ran;
   for (i = 0; i < 30000; i++) {
      h2p->Fill(50 * ran.Gaus(2., 1), ran.Gaus(2., 1), weight);
   }

   return h2p;
}

// test TH2Poly adding two histograms
TEST(TH2Poly, Add)
{
   // Create first hist
   TH2Poly *h2p_1 = createPoly(1);

   // Create second hist
   TH2Poly *h2p_2 = createPoly(0.1);

   // Create an added hist
   TH2Poly *h2p_added = (TH2Poly *)(h2p_1->Clone());
   h2p_added->Add(h2p_2, 1);

   EXPECT_LE(400, h2p_1->GetMaximum());
   EXPECT_LE(40, h2p_2->GetMaximum());
   EXPECT_LE(440, h2p_added->GetMaximum());

   delete h2p_1;
   delete h2p_2;
}

TEST(TH2Poly, CopyConstructor)
{
   TH2Poly *original = createPoly(1);
   TH2Poly copy{*original};

   original->Print();
   copy.Print();

   delete original;
}

TH2Poly *CreateHist()
{

   TH2Poly *h2p = new TH2Poly();
   Double_t x1[] = {0, 5, 6};
   Double_t y1[] = {0, 0, 5};
   Double_t x2[] = {0, -1, -1, 0};
   Double_t y2[] = {0, 0, -1, 3};
   Double_t x3[] = {4, 3, 0, 1, 2.4};
   Double_t y3[] = {4, 3.7, 1, 3.7, 2.5};
   h2p->AddBin(3, x1, y1);
   h2p->AddBin(4, x2, y2);
   h2p->AddBin(5, x3, y3);

   return h2p;
}

// test TH2Poly setting and retrieving bin error
TEST(TH2Poly, BinErrorUnweighted)
{

   auto h2p = CreateHist();
   // fill bin1
   for (int i = 0; i < 9; ++i)
      h2p->Fill(0.1, 0.01);
   EXPECT_EQ(3, h2p->GetBinError(1));

   // fill sea bins
   h2p->Fill(-0.5, -0.5, 1);
   h2p->Fill(2, -0.5, 1);
   EXPECT_EQ(sqrt(2), h2p->GetBinError(-5));

   // fill bin 3
   h2p->Fill(1, 3);
   h2p->Fill(1.5, 3);
   h2p->Fill(2.5, 3);
   EXPECT_EQ(sqrt(3), h2p->GetBinError(3));

   // fill overflow bin
   h2p->Fill(7, 3);
   EXPECT_EQ(1, h2p->GetBinError(-6));
   h2p->Fill(-2, 10);
   EXPECT_EQ(1, h2p->GetBinError(-1));

   h2p->SetBinContent(1, 10);
   EXPECT_EQ(sqrt(10), h2p->GetBinError(1));

   h2p->SetBinContent(-3, 4);
   EXPECT_EQ(2, h2p->GetBinError(-3));

   EXPECT_EQ(0, h2p->GetBinError(2));
   EXPECT_EQ(0, h2p->GetBinError(-2));
   EXPECT_EQ(0, h2p->GetBinError(-4));
   EXPECT_EQ(0, h2p->GetBinError(-7));
   EXPECT_EQ(0, h2p->GetBinError(-8));
   EXPECT_EQ(0, h2p->GetBinError(-9));
}

TEST(TH2Poly, BinErrorWeighted)
{
   auto h2p = CreateHist();

   // fill bins
   double w2 = 0;
   for (int i = 0; i < 10; ++i) {
      double w = gRandom->Uniform(0, 10);
      h2p->Fill(0.1, 0.01, w);
      w2 += w * w;
   }
   EXPECT_EQ(sqrt(w2), h2p->GetBinError(1));

   // fill sea bins
   h2p->Fill(-0.5, -0.5, 4);
   h2p->Fill(2, -0.5, 3);
   EXPECT_EQ(5, h2p->GetBinError(-5));

   // fill bin 3
   h2p->Fill(1, 3, 2);
   h2p->Fill(1.5, 3, 2);
   h2p->Fill(2.5, 3, 1);
   EXPECT_EQ(3, h2p->GetBinError(3));

   // fill overflow bin
   h2p->Fill(7, 3, 2);
   EXPECT_EQ(2, h2p->GetBinError(-6));
   h2p->Fill(-2, 10, 3);
   EXPECT_EQ(3, h2p->GetBinError(-1));

   EXPECT_EQ(0, h2p->GetBinError(2));
   EXPECT_EQ(0, h2p->GetBinError(-2));
   EXPECT_EQ(0, h2p->GetBinError(-3));
   EXPECT_EQ(0, h2p->GetBinError(-4));
   EXPECT_EQ(0, h2p->GetBinError(-7));
   EXPECT_EQ(0, h2p->GetBinError(-8));
   EXPECT_EQ(0, h2p->GetBinError(-9));
}

TEST(TH2, SetBinError)
{
   auto h2p = CreateHist();

   h2p->SetBinContent(1, 1.0);
   h2p->SetBinContent(2, 2.0);
   h2p->SetBinContent(3, 3.0);
   h2p->SetBinContent(-9, 4.0);

   EXPECT_EQ(1, h2p->GetBinError(1));
   EXPECT_EQ(sqrt(2), h2p->GetBinError(2));
   EXPECT_EQ(sqrt(3), h2p->GetBinError(3));
   EXPECT_EQ(2, h2p->GetBinError(-9));

   h2p->SetBinError(1, 1.5);
   h2p->SetBinError(2, 2.5);
   h2p->SetBinError(3, 3.5);
   h2p->SetBinError(-8, 1);

   EXPECT_EQ(1.5, h2p->GetBinError(1));
   EXPECT_EQ(2.5, h2p->GetBinError(2));
   EXPECT_EQ(3.5, h2p->GetBinError(3));
   EXPECT_EQ(1, h2p->GetBinError(-8));

   // setting a new content does not set bin error
   h2p->SetBinContent(-1, 3);
   EXPECT_EQ(0, h2p->GetBinError(-1));
}
