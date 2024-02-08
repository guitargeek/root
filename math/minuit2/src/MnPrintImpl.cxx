#include "MyMinuit2/MnPrint.h"

using ROOT::MyMinuit2::MnPrint;

#ifndef USE_ROOT_ERROR

#include <iostream>

#ifndef MN_OS
#define MN_OS std::cerr
#endif

void MnPrint::Impl(MnPrint::Verbosity level, const std::string &s)
{
   const char *label[4] = {"[Error]", "[Warn]", "[Info]", "[Debug]"};
   const int ilevel = static_cast<int>(level);
   MN_OS << label[ilevel] << " " << s << std::endl;
}

#else // use ROOT error reporting system

#include "TError.h"
#include <sstream>

void MnPrint::Impl(MnPrint::Verbosity level, const std::string &s)
{
   switch (level) {
   case MnPrint::eError: ::Error("MyMinuit2", "%s", s.c_str()); break;
   case MnPrint::eWarn: ::Warning("MyMinuit2", "%s", s.c_str()); break;
   case MnPrint::eInfo:  ::Info("MyMinuit2", "%s", s.c_str()); break;
   case MnPrint::eDebug: ::Info("MyMinuit2", "%s", s.c_str()); break;
   case MnPrint::eTrace: ::Info("MyMinuit2", "%s", s.c_str()); break;
   }
}

#endif // USE_ROOT_ERROR
