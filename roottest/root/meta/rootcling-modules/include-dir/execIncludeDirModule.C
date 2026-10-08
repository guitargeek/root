int execIncludeDirModule()
{
   TClass *cl = TClass::GetClass("IncludeDirModule");
   if (!cl || !cl->HasDictionary()) {
      std::cout << "No dictionary for IncludeDirModule\n";
      return 1;
   }
   std::cout << "fI = " << gROOT->ProcessLine("IncludeDirModule().fI;") << '\n';
   return 0;
}
