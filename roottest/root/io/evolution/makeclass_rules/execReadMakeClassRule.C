// Read back a TTree written with version 1 of MakeClassRule using version 2,
// where I/O customization rules apply, accessing the sub-branches directly in
// MakeClass mode (see https://github.com/root-project/root/issues/22611).

#include <TBranch.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdio>
#include <memory>

namespace {

// Surround the value with guards to detect writes to the wrong address.
struct Guarded {
   Int_t fBefore = -1;
   Int_t fValue = -1;
   Int_t fAfter = -1;
};

int check(const char *what, Long64_t entry, const Guarded &g, Int_t expected)
{
   std::printf("%s entry %lld: %d\n", what, entry, g.fValue);
   if (g.fValue != expected || g.fBefore != -1 || g.fAfter != -1) {
      std::printf("ERROR: expected %d, got %d (guards: %d, %d)\n", expected, g.fValue, g.fBefore, g.fAfter);
      return 1;
   }
   return 0;
}

} // namespace

int execReadMakeClassRule()
{
   if (gSystem->Load("libRoottestIoEvolutionMakeClassRuleV2") < 0)
      return 1;

   auto file = std::unique_ptr<TFile>(TFile::Open("makeclass_rules.root"));
   int errors = 0;

   // Branch that is the target of a rule, read on its own.
   {
      auto tree = file->Get<TTree>("t");
      Guarded run;
      TBranch *br = nullptr;
      tree->SetBranchAddress("obj.fRun", &run.fValue, &br);
      for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
         br->GetEntry(i);
         errors += check("fRun (branch)", i, run, 100 * (i + 1) + 1000);
      }
      delete tree;
   }

   // Branch that is the source of a rule, read on its own.
   {
      auto tree = file->Get<TTree>("t");
      Guarded src;
      TBranch *br = nullptr;
      tree->SetBranchAddress("obj.fSrc", &src.fValue, &br);
      for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
         br->GetEntry(i);
         errors += check("fSrc (branch)", i, src, 10 * (i + 1));
      }
      delete tree;
   }

   // Both branches, read via the tree.
   {
      auto tree = file->Get<TTree>("t");
      tree->SetMakeClass(1);
      Guarded run;
      Guarded src;
      tree->SetBranchAddress("obj.fRun", &run.fValue);
      tree->SetBranchAddress("obj.fSrc", &src.fValue);
      for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
         tree->GetEntry(i);
         errors += check("fRun (tree)", i, run, 100 * (i + 1) + 1000);
         errors += check("fSrc (tree)", i, src, 10 * (i + 1));
      }
      delete tree;
   }

   return errors;
}
