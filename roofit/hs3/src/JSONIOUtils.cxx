#include "JSONIOUtils.h"

#include <cmath>

bool startsWith(std::string_view str, std::string_view prefix)
{
   return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}

bool endsWith(std::string_view str, std::string_view suffix)
{
   return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

double roundPrecision(double d, int nSig)
{
   if (d == 0.0)
      return 0.0;
   int ndigits = std::floor(std::log10(std::abs(d))) + 1 - nSig;
   double sf = std::pow(10, ndigits);
   if (std::abs(d / sf) < 2)
      ndigits--;
   return sf * std::round(d / sf);
}
