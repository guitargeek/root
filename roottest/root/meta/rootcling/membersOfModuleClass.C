#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Load the dictionary `dict` and check that it provides the dictionaries for
// the comma-separated list of classes in `classNames`, and that the classes in
// `autoStreamerClassNames` use the automatic streamer, i.e. that they are
// selected by a rule with "+" although they have a ClassDef.
int membersOfModuleClass(const char *dict, const char *classNames, const char *autoStreamerClassNames = "")
{
   if (gSystem->Load(dict) < 0) {
      std::cerr << "Failed to load " << std::quoted(dict) << std::endl;
      return 1;
   }

   std::istringstream classNamesStream{classNames};
   std::string clName;
   while (std::getline(classNamesStream, clName, ',')) {
      auto cl = TClass::GetClass(clName.c_str());
      if (!cl) {
         std::cerr << "Class " << std::quoted(clName) << " not found!" << std::endl;
         return 2;
      }
      if (!cl->HasDictionary()) {
         std::cerr << "Class " << std::quoted(clName) << " has no dictionary!" << std::endl;
         return 3;
      }
   }

   std::istringstream autoStreamerClassNamesStream{autoStreamerClassNames};
   while (std::getline(autoStreamerClassNamesStream, clName, ',')) {
      auto cl = TClass::GetClass(clName.c_str());
      if (!cl || !cl->HasDictionary()) {
         std::cerr << "Class " << std::quoted(clName) << " not found or has no dictionary!" << std::endl;
         return 4;
      }
      if (cl->HasCustomStreamerMember()) {
         std::cerr << "Class " << std::quoted(clName) << " doesn't use the automatic streamer!" << std::endl;
         return 5;
      }
   }
   return 0;
}
