// Tests for the RooCmdArg
// Author: Jonas Rembser, CERN 11/2026

#include <RooArgSet.h>
#include <RooCmdArg.h>
#include <RooGlobalFunc.h>
#include <RooRealVar.h>

#include <gtest/gtest.h>

// Assigning a command argument gives the same as copy-constructing it.
TEST(RooCmdArg, Assignment)
{
   RooRealVar x("x", "x", 0, 10);
   RooRealVar y("y", "y", 0, 10);

   // payload sets are only there if the assigned command argument has them,
   // e.g. VisualizeError() without a parameter set must not get an empty one
   RooCmdArg arg = RooFit::Project(RooArgSet{x});
   ASSERT_NE(arg.getSet(0), nullptr);
   arg = RooFit::Range("r");
   EXPECT_EQ(arg.getSet(0), nullptr);
   RooCmdArg empty;
   empty = RooFit::Range("r");
   EXPECT_EQ(empty.getSet(0), nullptr);

   // sub-arguments are replaced, not appended
   RooCmdArg yVar = RooFit::YVar(y, RooFit::Binning(20));
   yVar = RooFit::YVar(y, RooFit::Binning(20));
   EXPECT_EQ(yVar.subArgs().GetSize(), 1);
}
