## \file
## \ingroup tutorial_roofit
## \notebook
##
## 'DATA AND CATEGORIES' RooFit tutorial macro #401
##
## Overview of advanced option for importing data from ROOT ROOT.TTree and ROOT.THx histograms
## Basic import options are demonstrated in rf102_dataimport.py
##
## \macro_code
## \macro_output
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C version)

import ROOT


# Import multiple TH1 into a RooDataHist
# ----------------------------------------------------------

# Create thee ROOT ROOT.TH1 histograms
hh_1 = ROOT.TH1D("hh1", "hh1", 100, -10., 10., 0, 3)
hh_2 = ROOT.TH1D("hh2", "hh2", 100, -10., 10., -3, 1)
hh_3 = ROOT.TH1D("hh3", "hh3", 100, -10., 10., +3, 4)

for i in range(1000):
    hh_1.Fill(ROOT.gRandom.Gaus(0., 3.))
    hh_2.Fill(ROOT.gRandom.Gaus(-3., 1.))
    hh_3.Fill(ROOT.gRandom.Gaus(+3., 4.))

# Declare observable x
x = ROOT.RooRealVar("x", "x", -10, 10)

# Create category observable c that serves as index for the ROOT histograms
c = ROOT.RooCategory("c", "c")
c.defineType("SampleA")
c.defineType("SampleB")
c.defineType("SampleC")
c.Print()

# Create a binned dataset that imports contents of all ROOT.TH1 mapped by
# index category c
dh = ROOT.RooDataHist("dh", "dh", [x], Index=c, Import={"SampleA": hh_1, "SampleB": hh_2, "SampleC": hh_3})
dh.Print()

dh2 = ROOT.RooDataHist("dh", "dh", [x], Index=c, Import={"SampleA": hh_1, "SampleB": hh_2, "SampleC": hh_3})
dh2.Print()

# Importing a ROOT TTree into a RooDataSet with cuts
# --------------------------------------------------------------------

treename = "tree"
filename = "rf401_importttreethx_py.root"

ROOT.RDataFrame(100) \
    .Define("x", "gRandom->Gaus(0, 3)") \
    .Define("y", "gRandom->Uniform() * 30 - 15") \
    .Define("z", "gRandom->Gaus(0, 5)") \
    .Define("i", "rdfentry_ % 3") \
    .Snapshot(treename, filename)

file = ROOT.TFile.Open(filename)
tree = file[treename]

# Define observables y,z
y = ROOT.RooRealVar("y", "y", -10, 10)
z = ROOT.RooRealVar("z", "z", -10, 10)

# Import only observables (y,z)
ds = ROOT.RooDataSet("ds", "ds", {x, y}, Import=tree)
ds.Print()

# Import observables (x,y,z) but only event for which (y+z<0) is ROOT.True
# Import observables (x,y,z) but only event for which (y+z<0) is ROOT.True
ds2 = ROOT.RooDataSet("ds2", "ds2", {x, y, z}, Import=tree, Cut="y+z<0")
ds2.Print()

# Importing integer ROOT TTree branches
# ---------------------------------------------------------------

# Import integer tree branch as ROOT.RooRealVar
i = ROOT.RooRealVar("i", "i", 0, 5)
ds3 = ROOT.RooDataSet("ds3", "ds3", {i, x}, Import=tree)
ds3.Print()

# Define category i
icat = ROOT.RooCategory("i", "i", {"State0": 0, "State1": 1})

# Import integer tree branch as ROOT.RooCategory (only events with i==0 and i==1
# will be imported as those are the only defined states)
ds4 = ROOT.RooDataSet("ds4", "ds4", {icat, x}, Import=tree)
ds4.Print()

# No need for the TTree anymore, so we can close the file that contains it
file.Close()

# Import multiple RooDataSets into a RooDataSet
# ----------------------------------------------------------------------------------------

# Create three ROOT.RooDataSets in (y,z)
dsA = ds2.reduce({x, y}, "z<-5")
dsB = ds2.reduce({x, y}, "abs(z)<5")
dsC = ds2.reduce({x, y}, "z>5")

# Create a dataset that imports contents of all the above datasets mapped
# by index category c
dsABC = ROOT.RooDataSet("dsABC", "dsABC", {x, y}, Index=c, Import={"SampleA": dsA, "SampleB": dsB, "SampleC": dsC})

dsABC.Print()
