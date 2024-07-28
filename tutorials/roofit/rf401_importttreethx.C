/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Data and categories: advanced options for importing data from ROOT TTree and THx histograms
///
/// Basic import options are demonstrated in rf102_dataimport.C
///
/// \macro_code
/// \macro_output
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "ROOT/RDataFrame.hxx"

#include "TAxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"

#include <map>

TH1 *makeTH1(const char *name, double mean, double sigma);

void rf401_importttreethx()
{
   using namespace RooFit;

   // I m p o r t  m u l t i p l e   T H 1   i n t o   a   R o o D a t a H i s t
   // --------------------------------------------------------------------------

   // Create thee ROOT TH1 histograms
   TH1 *hh_1 = makeTH1("hh1", 0, 3);
   TH1 *hh_2 = makeTH1("hh2", -3, 1);
   TH1 *hh_3 = makeTH1("hh3", +3, 4);

   // Declare observable x
   RooRealVar x("x", "x", -10, 10);

   // Create category observable c that serves as index for the ROOT histograms
   RooCategory c("c", "c", {{"SampleA",0}, {"SampleB",1}, {"SampleC",2}});

   // Create a binned dataset that imports contents of all TH1 mapped by index category c
   RooDataHist *dh = new RooDataHist("dh", "dh", x, Index(c), Import("SampleA", *hh_1), Import("SampleB", *hh_2),
                                     Import("SampleC", *hh_3));
   dh->Print();

   // Alternative constructor form for importing multiple histograms
   std::map<std::string, TH1 *> hmap;
   hmap["SampleA"] = hh_1;
   hmap["SampleB"] = hh_2;
   hmap["SampleC"] = hh_3;
   RooDataHist *dh2 = new RooDataHist("dh", "dh", x, c, hmap);
   dh2->Print();

   // I m p o r t i n g   a   T T r e e   i n t o   a   R o o D a t a S e t   w i t h   c u t s
   // -----------------------------------------------------------------------------------------

   const std::string treename = "tree";
   const std::string filename = "rf401_importttreethx.root";

   ROOT::RDataFrame(100)
       .Define("x", "gRandom->Gaus(0, 3)")
       .Define("y", "gRandom->Uniform() * 30 - 15")
       .Define("z", "gRandom->Gaus(0, 5)")
       .Define("i", "rdfentry_ % 3")
       .Snapshot(treename, filename);

   std::unique_ptr<TFile> file{TFile::Open(filename.c_str())};
   TTree *tree = file->Get<TTree>(treename.c_str());

   // Define observables y,z
   RooRealVar y("y", "y", -10, 10);
   RooRealVar z("z", "z", -10, 10);

   // Import only observables (y,z)
   RooDataSet ds("ds", "ds", {x, y}, Import(*tree));
   ds.Print();

   // Import observables (x,y,z) but only event for which (y+z<0) is true
   RooDataSet ds2("ds2", "ds2", {x, y, z}, Import(*tree), Cut("y+z<0"));
   ds2.Print();

   // I m p o r t i n g   i n t e g e r   T T r e e   b r a n c h e s
   // ---------------------------------------------------------------

   // Import integer tree branch as RooRealVar
   RooRealVar i("i", "i", 0, 5);
   RooDataSet ds3("ds3", "ds3", {i, x}, Import(*tree));
   ds3.Print();

   // Define category i
   RooCategory icat("i", "i");
   icat.defineType("State0", 0);
   icat.defineType("State1", 1);

   // Import integer tree branch as RooCategory (only events with i==0 and i==1
   // will be imported as those are the only defined states)
   RooDataSet ds4("ds4", "ds4", {icat, x}, Import(*tree));
   ds4.Print();

   // No need for the TTree anymore, so we can close the file that contains it
   file->Close();

   // I m p o r t  m u l t i p l e   R o o D a t a S e t s   i n t o   a   R o o D a t a S e t
   // ----------------------------------------------------------------------------------------

   // Create three RooDataSets in (y,z)
   std::unique_ptr<RooAbsData> dsA{ds2.reduce({x, y}, "z<-5")};
   std::unique_ptr<RooAbsData> dsB{ds2.reduce({x, y}, "abs(z)<5")};
   std::unique_ptr<RooAbsData> dsC{ds2.reduce({x, y}, "z>5")};

   // Create a dataset that imports contents of all the above datasets mapped by index category c
   RooDataSet dsABC{"dsABC", "dsABC", {x, y}, Index(c), Import("SampleA", *dsA),
                    Import("SampleB", *dsB), Import("SampleC", *dsC)};

   dsABC.Print();
}

TH1 *makeTH1(const char *name, double mean, double sigma)
{
   // Create ROOT TH1 filled with a Gaussian distribution

   TH1D *hh = new TH1D(name, name, 100, -10, 10);
   for (int i = 0; i < 1000; i++) {
      hh->Fill(gRandom->Gaus(mean, sigma));
   }
   return hh;
}
