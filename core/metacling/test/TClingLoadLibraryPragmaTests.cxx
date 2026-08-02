/// \file TClingLoadLibraryPragmaTests.cxx
/// Regression tests for resolving libraries passed to `#pragma cling load`,
/// which is what the R__LOAD_LIBRARY macro expands to. In particular, an
/// absolute path with the shared-library extension omitted must resolve, just
/// like a relative path does. See DynamicLibraryManager::lookupLibrary.

#include "TInterpreter.h"
#include "TSystem.h"

#include <ROOT/TestSupport.hxx>

#include "gtest/gtest.h"

#include <fstream>
#include <string>

namespace {

// Compile a tiny shared library exporting `extern "C" int <symbol>()` (returning
// 42) WITHOUT loading it (ACLiC option "c"), and return the absolute path to the
// produced library with the shared-library extension stripped off.
std::string CompileUnloadedHelper(const std::string &symbol)
{
   const std::string dir = gSystem->TempDirectory();
   const std::string src = dir + "/" + symbol + ".cxx";
   {
      std::ofstream srcFile(src);
      srcFile << "extern \"C\" int " << symbol << "() { return 42; }\n";
   }

   // "k" keep the library, "f" force a fresh build, "c" compile but do not load.
   EXPECT_EQ(gSystem->CompileMacro(src.c_str(), "kfc"), 1) << "Failed to compile " << src;

   // ACLiC names the library "<stem>_cxx.<soext>" next to the source file.
   const std::string stem = dir + "/" + symbol + "_cxx";
   EXPECT_FALSE(gSystem->AccessPathName((stem + "." + gSystem->GetSoExt()).c_str()))
      << "Expected library was not produced for " << src;
   return stem;
}

long CallSymbol(const std::string &symbol, TInterpreter::EErrorCode &err)
{
   err = TInterpreter::kNoError;
   gInterpreter->Declare(("extern \"C\" int " + symbol + "();").c_str());
   return gInterpreter->ProcessLine((symbol + "();").c_str(), &err);
}

} // namespace

// An absolute path with the shared-library extension omitted must resolve.
// Before the DynamicLibraryManager fix this failed with a "file not found"
// error, because absolute paths skipped the extension-adding lookup that
// relative paths went through.
TEST(TClingLoadLibraryPragma, AbsolutePathWithoutExtension)
{
   const std::string symbol = "cling_load_abs_noext";
   const std::string stem = CompileUnloadedHelper(symbol);

   gInterpreter->ProcessLine(("#pragma cling load(\"" + stem + "\")").c_str());

   TInterpreter::EErrorCode err;
   EXPECT_EQ(CallSymbol(symbol, err), 42);
   EXPECT_EQ(err, TInterpreter::kNoError);
}

// Baseline: the same absolute path WITH the extension has always worked.
TEST(TClingLoadLibraryPragma, AbsolutePathWithExtension)
{
   const std::string symbol = "cling_load_abs_ext";
   const std::string lib = CompileUnloadedHelper(symbol) + "." + gSystem->GetSoExt();

   gInterpreter->ProcessLine(("#pragma cling load(\"" + lib + "\")").c_str());

   TInterpreter::EErrorCode err;
   EXPECT_EQ(CallSymbol(symbol, err), 42);
   EXPECT_EQ(err, TInterpreter::kNoError);
}

// A nonexistent absolute path must still be rejected: the extension-adding
// lookup should not turn a missing library into a false positive. The pragma is
// expected to report "file not found".
TEST(TClingLoadLibraryPragma, NonexistentAbsolutePath)
{
   const std::string stem = std::string(gSystem->TempDirectory()) + "/cling_load_abs_missing";

   ROOT::TestSupport::CheckDiagsRAII diags;
   diags.requiredDiag(kError, "cling", "file not found", /*matchFullMessage=*/false);

   gInterpreter->ProcessLine(("#pragma cling load(\"" + stem + "\")").c_str());
}
