// Read data (CERN staff) from an ascii file and create a root file with a Tree.
// see also a variant in staff.C
// Author: Rene Brun
   
TFile *cernbuild(Int_t get=0) {

   Int_t           Category;
   UInt_t          Flag;
   Int_t           Age;
   Int_t           Service;
   Int_t           Children;
   Int_t           Grade;
   Int_t           Step;
   Int_t           Hrweek;
   Int_t           Cost;
   Char_t          Division[4];
   Char_t          Nation[3];

   //The input file cern.dat is a copy of the CERN staff data base
   //from 1988
   TString filename = "cernstaff.root";
   FILE *fp = fopen("cernstaff.dat","r");

   TFile *hfile = 0;
   if (get) {
      // if the argument get =1 return the file "cernstaff.root"
      // if the file does not exist, it is created
      if (!gSystem->AccessPathName("cernstaff.root",kFileExists)) {
         hfile = TFile::Open("cernstaff.root"); //in $ROOTSYS/tutorials/io/tree
         if (hfile) return hfile;
      }
   }
   //no cernstaff.root file found. Must generate it !
   //generate cernstaff.root in $ROOTSYS/tutorials/io/tree if we have write access
   if (!gSystem->AccessPathName(".",kWritePermission)) {
      //otherwise generate cernstaff.root in the current directory
   } else {
      printf("you must run the script in a directory with write access\n");
      return 0;
   }
   hfile = TFile::Open(filename,"RECREATE");
   TTree *tree = new TTree("T","CERN 1988 staff data");
   tree->Branch("Category",&Category,"Category/I");
   tree->Branch("Flag",&Flag,"Flag/i");
   tree->Branch("Age",&Age,"Age/I");
   tree->Branch("Service",&Service,"Service/I");
   tree->Branch("Children",&Children,"Children/I");
   tree->Branch("Grade",&Grade,"Grade/I");
   tree->Branch("Step",&Step,"Step/I");
   tree->Branch("Hrweek",&Hrweek,"Hrweek/I");
   tree->Branch("Cost",&Cost,"Cost/I");
   tree->Branch("Division",Division,"Division/C");
   tree->Branch("Nation",Nation,"Nation/C");
   char line[80];
   while (fgets(line,80,fp)) {
      sscanf(&line[0],"%d %d %d %d %d",&Category,&Flag,&Age,&Service,&Children);
      sscanf(&line[32],"%d %d  %d %d %s %s",&Grade,&Step,&Hrweek,&Cost,Division,Nation);
      tree->Fill();
   }
   // tree->Print();
   tree->Write();

   fclose(fp);
   delete hfile;
   if (get) {
      //we come here when the script is executed outside $ROOTSYS/tutorials/io/tree
      hfile = TFile::Open(filename);
      return hfile;
   }
   return 0;
}

int execDrawString()
{
   cernbuild(1);
#ifdef ClingWorkAroundMissingDynamicScope
   TTree *T; gFile->GetObject("T",T);
#endif
   T->Draw("Division");
#ifdef ClingWorkAroundMissingDynamicScope
   TH1F *htemp = (TH1F*)gROOT->FindObject("htemp");
#endif
   htemp->Print("all");
#ifdef ClingWorkAroundErracticValuePrinter
   printf("(int)0\n");
#endif
   return 0;
}

