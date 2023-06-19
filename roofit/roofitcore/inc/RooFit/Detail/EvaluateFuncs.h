/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2023
 *   Garima Singh, CERN 2023
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_Detail_EvaluateFuncs_h
#define RooFit_Detail_EvaluateFuncs_h

#include <TMath.h>

#include <cmath>

namespace RooFit {

namespace Detail {

namespace EvaluateFuncs {

/// @brief Function to evaluate an un-normalized RooGaussian.
inline double gaussianEvaluate(double x, double mean, double sigma)
{
   const double arg = x - mean;
   const double sig = sigma;
   return std::exp(-0.5 * arg * arg / (sig * sig));
}

/// In pdfMode, a coefficient for the constant term of 1.0 is implied if lowestOrder > 0.
template <bool pdfMode = false>
inline double polynomialEvaluate(double const *coeffs, int nCoeffs, int lowestOrder, double x)
{
   double retVal = coeffs[nCoeffs - 1];
   for (int i = nCoeffs - 2; i >= 0; i--)
      retVal = coeffs[i] + x * retVal;
   retVal = retVal * std::pow(x, lowestOrder);
   return retVal + (pdfMode && lowestOrder > 0 ? 1.0 : 0.0);
}

inline double chebychevEvaluate(double *coeffs, unsigned int nCoeffs, double x_in, double xMin, double xMax)
{
   // transform to range [-1, +1]
   const double xPrime = (x_in - 0.5 * (xMax + xMin)) / (0.5 * (xMax - xMin));

   // extract current values of coefficients
   double sum = 1.;
   if (nCoeffs > 0) {
      double curr = xPrime;
      double twox = 2 * xPrime;
      double last = 1;
      double newval = twox * curr - last;
      last = curr;
      curr = newval;
      for (unsigned int i = 0; nCoeffs != i; ++i) {
         sum += last * coeffs[i];
         newval = twox * curr - last;
         last = curr;
         curr = newval;
      }
   }
   return sum;
}

inline double constraintSumEvaluate(double const *comp, unsigned int compSize)
{
   double sum = 0;
   for (unsigned int i = 0; i < compSize; i++) {
      sum -= std::log(comp[i]);
   }
   return sum;
}

inline unsigned int getUniformBinning(double low, double high, double val, unsigned int numBins)
{
   double binWidth = (high - low) / numBins;
   return val >= high ? numBins - 1 : std::abs((val - low) / binWidth);
}

inline double poissonEvaluate(double x, double par)
{
   if (par < 0)
      return TMath::QuietNaN();

   if (x < 0)
      return 0;
   else if (x == 0.0)
      return std::exp(-par);
   else {
      double out = x * std::log(par) - TMath::LnGamma(x + 1.) - par;
      return std::exp(out);
   }
}

inline double polyInterpValue(unsigned int i, double *polCoeff, double const *low, double const *high, double x,
                              double boundary, double n, double nominal, unsigned int &logInit)
{
   // code for polynomial interpolation used when interpCode=4
   double x0 = boundary;

   // cache the polynomial coefficient values
   // which do not depend on x but on the boundaries values
   if (!logInit) {
      logInit = 1;
      for (unsigned int j = 0; j < n; j++) {
         // location of the 6 coefficient for the j-th variable
         unsigned int offset = j * 6;

         // GHL: Swagato's suggestions
         double pow_up = std::pow(high[j] / nominal, x0);
         double pow_down = std::pow(low[j] / nominal, x0);
         double logHi = std::log(high[j]);
         double logLo = std::log(low[j]);
         double pow_up_log = high[j] <= 0.0 ? 0.0 : pow_up * logHi;
         double pow_down_log = low[j] <= 0.0 ? 0.0 : -pow_down * logLo;
         double pow_up_log2 = high[j] <= 0.0 ? 0.0 : pow_up_log * logHi;
         double pow_down_log2 = low[j] <= 0.0 ? 0.0 : -pow_down_log * logLo;

         double S0 = (pow_up + pow_down) / 2;
         double A0 = (pow_up - pow_down) / 2;
         double S1 = (pow_up_log + pow_down_log) / 2;
         double A1 = (pow_up_log - pow_down_log) / 2;
         double S2 = (pow_up_log2 + pow_down_log2) / 2;
         double A2 = (pow_up_log2 - pow_down_log2) / 2;

         // cache  coefficient of the polynomial
         polCoeff[0 + offset] = 1. / (8 * x0) * (15 * A0 - 7 * x0 * S1 + x0 * x0 * A2);
         polCoeff[1 + offset] = 1. / (8 * x0 * x0) * (-24 + 24 * S0 - 9 * x0 * A1 + x0 * x0 * S2);
         polCoeff[2 + offset] = 1. / (4 * std::pow(x0, 3)) * (-5 * A0 + 5 * x0 * S1 - x0 * x0 * A2);
         polCoeff[3 + offset] = 1. / (4 * std::pow(x0, 4)) * (12 - 12 * S0 + 7 * x0 * A1 - x0 * x0 * S2);
         polCoeff[4 + offset] = 1. / (8 * std::pow(x0, 5)) * (+3 * A0 - 3 * x0 * S1 + x0 * x0 * A2);
         polCoeff[5 + offset] = 1. / (8 * std::pow(x0, 6)) * (-8 + 8 * S0 - 5 * x0 * A1 + x0 * x0 * S2);
      }
   }

   unsigned int offset = 6 * i;

   double a = polCoeff[offset + 0];
   double b = polCoeff[offset + 1];
   double c = polCoeff[offset + 2];
   double d = polCoeff[offset + 3];
   double e = polCoeff[offset + 4];
   double f = polCoeff[offset + 5];

   // evaluate the 6-th degree polynomial using Horner's method
   return 1. + x * (a + x * (b + x * (c + x * (d + x * (e + x * f)))));
}

inline double flexibleInterpVarProcessParam(unsigned int code, unsigned int i, double *polCoeff, double const *low,
                                            double const *high, double n, double boundary, double nominal,
                                            double paramVal, double total, unsigned int &logInit)
{
   if (code == 0) {
      // piece-wise linear
      if (paramVal > 0)
         return total + paramVal * (high[i] - nominal);
      else
         return total + paramVal * (nominal - low[i]);
   } else if (code == 1) {
      // pice-wise log
      if (paramVal >= 0)
         return total * pow(high[i] / nominal, +paramVal);
      else
         return total * pow(low[i] / nominal, -paramVal);
   } else if (code == 2) {
      // parabolic with linear
      double a = 0.5 * (high[i] + low[i]) - nominal;
      double b = 0.5 * (high[i] - low[i]);
      double c = 0;
      if (paramVal > 1) {
         return total + (2 * a + b) * (paramVal - 1) + high[i] - nominal;
      } else if (paramVal < -1) {
         return total + -1 * (2 * a - b) * (paramVal + 1) + low[i] - nominal;
      } else {
         return total + a * pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 3) {
      // parabolic version of log-normal
      double a = 0.5 * (high[i] + low[i]) - nominal;
      double b = 0.5 * (high[i] - low[i]);
      double c = 0;
      if (paramVal > 1) {
         return total + (2 * a + b) * (paramVal - 1) + high[i] - nominal;
      } else if (paramVal < -1) {
         return total + -1 * (2 * a - b) * (paramVal + 1) + low[i] - nominal;
      } else {
         return total + a * pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 4) {
      double x = paramVal;
      if (x >= boundary) {
         return total * std::pow(high[i] / nominal, +paramVal);
      } else if (x <= -boundary) {
         return total * std::pow(low[i] / nominal, -paramVal);
      } else if (x != 0) {
         return total * polyInterpValue(i, polCoeff, low, high, x, boundary, n, nominal, logInit);
      }
   }

   return 0;
}

inline double flexibleInterpVarEvaluate(unsigned int const *code, double const *params, unsigned int paramSize,
                                        double const *low, double const *high, double nominal, double boundary)
{
   double polCoeff[6 * paramSize];
   double total = nominal;
   unsigned int logInit = 0;
   for (std::size_t i = 0; i < paramSize; ++i) {
      total = flexibleInterpVarProcessParam(code[i], i, polCoeff, low, high, paramSize, boundary, nominal, params[i],
                                            total, logInit);
   }
   if (total <= 0) {
      total = TMath::Limits<double>::Min();
   }
   return total;
}

} // namespace EvaluateFuncs

} // namespace Detail

} // namespace RooFit

#endif
