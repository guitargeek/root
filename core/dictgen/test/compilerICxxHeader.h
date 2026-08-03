// Header used by the "rootcling-compilerI-path-fallback" regression test.
//
// It pulls in a C++ standard library header. If rootcling resolves the C++
// standard library against the wrong compiler (e.g. a stale one found on $PATH
// instead of the one driving the build, passed via -compilerI), parsing this
// header fails -- mirroring the original report where <concepts> could not be
// found. See the test driver compilerIPathFallback.cmake for details.

#ifndef ROOTTEST_DICTGEN_COMPILERI_CXXHEADER_H
#define ROOTTEST_DICTGEN_COMPILERI_CXXHEADER_H

#include <vector>

class CompilerIProbe {
public:
   std::vector<int> fValues;
};

#endif
