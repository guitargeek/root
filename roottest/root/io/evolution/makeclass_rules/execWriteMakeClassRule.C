// Write a TTree with version 1 of MakeClassRule, to be read back with
// version 2 by execReadMakeClassRule.C.

#include <TClass.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <memory>

int execWriteMakeClassRule()
{
   if (gSystem->Load("libRoottestIoEvolutionMakeClassRuleV1") < 0)
      return 1;

   TClass *cl = TClass::GetClass("MakeClassRule");
   void *obj = cl->New();
   auto setMember = [&](const char *name, Int_t value) {
      *reinterpret_cast<Int_t *>(static_cast<char *>(obj) + cl->GetDataMemberOffset(name)) = value;
   };

   auto file = std::unique_ptr<TFile>(TFile::Open("makeclass_rules.root", "RECREATE"));
   TTree tree("t", "t");
   tree.Branch("obj", "MakeClassRule", &obj);
   for (Int_t i = 1; i <= 2; ++i) {
      setMember("fDropped", -100 * i);
      setMember("fRun", 100 * i);
      setMember("fSrc", 10 * i);
      tree.Fill();
   }
   tree.Write();
   file->Close();
   cl->Destructor(obj);
   return 0;
}
