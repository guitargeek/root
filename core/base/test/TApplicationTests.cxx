#include "TApplication.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

// https://github.com/root-project/root/issues/23648
TEST(TApplication, ExpressionsKeepCommandLineOrder)
{
   // Positional file arguments are only considered if the file exists and is not empty.
   const char *file1 = "TApplicationTests_1.root";
   const char *file2 = "TApplicationTests_2.root";
   for (const char *fname : {file1, file2})
      std::ofstream(fname) << "dummy";

   std::vector<std::string> args{
      "TApplicationTests", "-b", "-e", "first", file1, "--execute", "second", "-e", "third", file2, "-e", "fourth"};
   std::vector<char *> argv;
   for (auto &arg : args)
      argv.push_back(arg.data());
   argv.push_back(nullptr);
   int argc = args.size();

   TApplication app("TApplicationTests", &argc, argv.data());

   std::remove(file1);
   std::remove(file2);

   const std::vector<std::string> expected{"first", file1, "second", "third", file2, "fourth"};
   TObjArray *inputs = app.InputFiles();
   ASSERT_NE(inputs, nullptr);
   ASSERT_EQ(inputs->GetEntries(), static_cast<int>(expected.size()));
   for (std::size_t i = 0; i < expected.size(); ++i) {
      auto *input = static_cast<TObjString *>(inputs->At(i));
      const bool isExpression = input->TestBit(TApplication::kExpression);
      if (expected[i].find(".root") == std::string::npos) {
         EXPECT_TRUE(isExpression) << "input " << i;
         EXPECT_EQ(input->String(), expected[i]);
      } else {
         EXPECT_FALSE(isExpression) << "input " << i;
         EXPECT_TRUE(input->String().EndsWith(expected[i].c_str())) << input->String();
      }
   }
}
