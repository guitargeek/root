#include <iomanip>
#include <iostream>

int nestedTemplateOfModuleClass()
{
   gSystem->Load("nestedTemplateOfModuleClassDict");

   for (const char *clName :
        {"NestedTemplateArg", "ROOT::Detail::TCollectionProxyInfo::Environ<NestedTemplateArg>",
         "ROOT::Detail::TCollectionProxyInfo::PairHolder<NestedTemplateArg,int>"}) {
      auto cl = TClass::GetClass(clName);
      if (!cl) {
         std::cerr << "Class " << std::quoted(clName) << " not found!" << std::endl;
         return 1;
      }
      if (!cl->HasDictionary()) {
         std::cerr << "Class " << std::quoted(clName) << " has no dictionary!" << std::endl;
         return 2;
      }
   }
   return 0;
}
