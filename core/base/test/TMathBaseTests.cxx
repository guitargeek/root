#include "gtest/gtest.h"

#include "TMathBase.h"

#include <list>

// Tests for TMath::BinarySearch().
TEST(TMathBase, BinarySearch)
{
   // Use std::list for the test to make sure that the implementation doesn't
   // "cheat" by assuming contiguous memory.
   std::list<double> l1{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

   auto search = [&l1](double x) { return TMath::BinarySearch(l1.begin(), l1.end(), x); };

   // If the value is smaller than the smallest element in the list, return the end() iterator.
   EXPECT_EQ(search(-0.5), l1.end());

   for (double x : l1) {
      // Looking for the value should return the value.
      EXPECT_EQ(*search(x), x);
      // If no match found, function gives nearest element smaller than value.
      EXPECT_EQ(*search(x + 0.5), x);
   }
}
