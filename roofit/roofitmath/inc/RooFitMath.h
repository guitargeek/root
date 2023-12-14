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

#ifndef RooFitMath_h
#define RooFitMath_h

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>
#include <TMath.h>

#include <cmath>

namespace RooFitMath {

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

   if (x < 0) {
      return 0;
   } else if (x == 0.0) {
      return std::exp(-par);
   } else {
      double out = x * std::log(par) - TMath::LnGamma(x + 1.) - par;
      return std::exp(out);
   }
}

inline double interpolate6thDegree(double x, double low, double high, double nominal, double boundary)
{
   double t = x / boundary;
   double eps_plus = high - nominal;
   double eps_minus = nominal - low;
   double S = 0.5 * (eps_plus + eps_minus);
   double A = 0.0625 * (eps_plus - eps_minus);

   return x * (S + t * A * (15 + t * t * (-10 + t * t * 3)));
}

inline double interpolate6thDegreeExp(double x, double low, double high, double nominal, double boundary)
{
   double x0 = boundary;

   // GHL: Swagato's suggestions
   double powUp = std::pow(high / nominal, x0);
   double powDown = std::pow(low / nominal, x0);
   double logHi = std::log(high);
   double logLo = std::log(low);
   double powUpLog = high <= 0.0 ? 0.0 : powUp * logHi;
   double powDownLog = low <= 0.0 ? 0.0 : -powDown * logLo;
   double powUpLog2 = high <= 0.0 ? 0.0 : powUpLog * logHi;
   double powDownLog2 = low <= 0.0 ? 0.0 : -powDownLog * logLo;

   double S0 = 0.5 * (powUp + powDown);
   double A0 = 0.5 * (powUp - powDown);
   double S1 = 0.5 * (powUpLog + powDownLog);
   double A1 = 0.5 * (powUpLog - powDownLog);
   double S2 = 0.5 * (powUpLog2 + powDownLog2);
   double A2 = 0.5 * (powUpLog2 - powDownLog2);

   // fcns+der+2nd_der are eq at bd

   double a = 1. / (8 * x0) * (15 * A0 - 7 * x0 * S1 + x0 * x0 * A2);
   double b = 1. / (8 * x0 * x0) * (-24 + 24 * S0 - 9 * x0 * A1 + x0 * x0 * S2);
   double c = 1. / (4 * std::pow(x0, 3)) * (-5 * A0 + 5 * x0 * S1 - x0 * x0 * A2);
   double d = 1. / (4 * std::pow(x0, 4)) * (12 - 12 * S0 + 7 * x0 * A1 - x0 * x0 * S2);
   double e = 1. / (8 * std::pow(x0, 5)) * (+3 * A0 - 3 * x0 * S1 + x0 * x0 * A2);
   double f = 1. / (8 * std::pow(x0, 6)) * (-8 + 8 * S0 - 5 * x0 * A1 + x0 * x0 * S2);

   // evaluate the 6-th degree polynomial using Horner's method
   double value = 1. + x * (a + x * (b + x * (c + x * (d + x * (e + x * f)))));
   return value;
}

inline double
flexibleInterp(unsigned int code, double low, double high, double boundary, double nominal, double paramVal, double res)
{
   if (code == 0) {
      // piece-wise linear
      if (paramVal > 0) {
         return paramVal * (high - nominal);
      } else {
         return paramVal * (nominal - low);
      }
   } else if (code == 1) {
      // piece-wise log
      if (paramVal >= 0) {
         return res * (std::pow(high / nominal, +paramVal) - 1);
      } else {
         return res * (std::pow(low / nominal, -paramVal) - 1);
      }
   } else if (code == 2) {
      // parabolic with linear
      double a = 0.5 * (high + low) - nominal;
      double b = 0.5 * (high - low);
      double c = 0;
      if (paramVal > 1) {
         return (2 * a + b) * (paramVal - 1) + high - nominal;
      } else if (paramVal < -1) {
         return -1 * (2 * a - b) * (paramVal + 1) + low - nominal;
      } else {
         return a * std::pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 3) {
      // parabolic version of log-normal
      double a = 0.5 * (high + low) - nominal;
      double b = 0.5 * (high - low);
      double c = 0;
      if (paramVal > 1) {
         return (2 * a + b) * (paramVal - 1) + high - nominal;
      } else if (paramVal < -1) {
         return -1 * (2 * a - b) * (paramVal + 1) + low - nominal;
      } else {
         return a * std::pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 4) {
      double x = paramVal;
      if (x >= boundary) {
         return x * (high - nominal);
      } else if (x <= -boundary) {
         return x * (nominal - low);
      }

      return interpolate6thDegree(x, low, high, nominal, boundary);
   } else if (code == 5) {
      double x = paramVal;
      double mod = 1.0;
      if (x >= boundary) {
         mod = std::pow(high / nominal, +paramVal);
      } else if (x <= -boundary) {
         mod = std::pow(low / nominal, -paramVal);
      } else {
         mod = interpolate6thDegreeExp(x, low, high, nominal, boundary);
      }
      return res * (mod - 1.0);
   }

   return 0.0;
}

inline double logNormalEvaluate(double x, double k, double m0)
{
   return ROOT::Math::lognormal_pdf(x, std::log(m0), std::abs(std::log(k)));
}

inline double logNormalEvaluateStandard(double x, double sigma, double mu)
{
   return ROOT::Math::lognormal_pdf(x, mu, std::abs(sigma));
}

/// @brief Function to calculate the integral of an un-normalized RooGaussian over x. To calculate the integral over
/// mean, just interchange the respective values of x and mean.
/// @param xMin Minimum value of variable to integrate wrt.
/// @param xMax Maximum value of of variable to integrate wrt.
/// @param mean Mean.
/// @param sigma Sigma.
/// @return The integral of an un-normalized RooGaussian over the value in x.
inline double gaussianIntegral(double xMin, double xMax, double mean, double sigma)
{
   // The normalisation constant 1./sqrt(2*pi*sigma^2) is left out in evaluate().
   // Therefore, the integral is scaled up by that amount to make RooFit normalise
   // correctly.
   double resultScale = std::sqrt(TMath::TwoPi()) * sigma;

   // Here everything is scaled and shifted into a standard normal distribution:
   double xscale = TMath::Sqrt2() * sigma;
   double scaledMin = 0.;
   double scaledMax = 0.;
   scaledMin = (xMin - mean) / xscale;
   scaledMax = (xMax - mean) / xscale;

   // Here we go for maximum precision: We compute all integrals in the UPPER
   // tail of the Gaussian, because erfc has the highest precision there.
   // Therefore, the different cases for range limits in the negative hemisphere are mapped onto
   // the equivalent points in the upper hemisphere using erfc(-x) = 2. - erfc(x)
   double ecmin = TMath::Erfc(std::abs(scaledMin));
   double ecmax = TMath::Erfc(std::abs(scaledMax));

   double cond = 0.0;
   // Don't put this "prd" inside the "if" because clad will not be able to differentiate the code correctly (as of
   // v1.1)!
   double prd = scaledMin * scaledMax;
   if (prd < 0.0) {
      cond = 2.0 - (ecmin + ecmax);
   } else if (scaledMax <= 0.0) {
      cond = ecmax - ecmin;
   } else {
      cond = ecmin - ecmax;
   }
   return resultScale * 0.5 * cond;
}

inline double exponentialIntegral(double xMin, double xMax, double constant)
{
   if (constant == 0.0) {
      return xMax - xMin;
   }

   return (std::exp(constant * xMax) - std::exp(constant * xMin)) / constant;
}

/// In pdfMode, a coefficient for the constant term of 1.0 is implied if lowestOrder > 0.
template <bool pdfMode = false>
inline double polynomialIntegral(double const *coeffs, int nCoeffs, int lowestOrder, double xMin, double xMax)
{
   int denom = lowestOrder + nCoeffs;
   double min = coeffs[nCoeffs - 1] / double(denom);
   double max = coeffs[nCoeffs - 1] / double(denom);

   for (int i = nCoeffs - 2; i >= 0; i--) {
      denom--;
      min = (coeffs[i] / double(denom)) + xMin * min;
      max = (coeffs[i] / double(denom)) + xMax * max;
   }

   max = max * std::pow(xMax, 1 + lowestOrder);
   min = min * std::pow(xMin, 1 + lowestOrder);

   return max - min + (pdfMode && lowestOrder > 0.0 ? xMax - xMin : 0.0);
}

/// use fast FMA if available, fall back to normal arithmetic if not
inline double fast_fma(double x, double y, double z) noexcept
{
#if defined(FP_FAST_FMA) // check if std::fma has fast hardware implementation
   return std::fma(x, y, z);
#else // defined(FP_FAST_FMA)
   // std::fma might be slow, so use a more pedestrian implementation
#if defined(__clang__)
#pragma STDC FP_CONTRACT ON // hint clang that using an FMA is okay here
#endif                      // defined(__clang__)
   return (x * y) + z;
#endif                      // defined(FP_FAST_FMA)
}

inline double chebychevIntegral(double const *coeffs, unsigned int nCoeffs, double xMin, double xMax, double xMinFull,
                                double xMaxFull)
{
   const double halfrange = .5 * (xMax - xMin);
   const double mid = .5 * (xMax + xMin);

   // the full range of the function is mapped to the normalised [-1, 1] range
   const double b = (xMaxFull - mid) / halfrange;
   const double a = (xMinFull - mid) / halfrange;

   // coefficient for integral(T_0(x)) is 1 (implicit), integrate by hand
   // T_0(x) and T_1(x), and use for n > 1: integral(T_n(x) dx) =
   // (T_n+1(x) / (n + 1) - T_n-1(x) / (n - 1)) / 2
   double sum = b - a; // integrate T_0(x) by hand

   const unsigned int iend = nCoeffs;
   if (iend > 0) {
      {
         // integrate T_1(x) by hand...
         const double c = coeffs[0];
         sum = fast_fma(0.5 * (b + a) * (b - a), c, sum);
      }
      if (1 < iend) {
         double bcurr = b;
         double btwox = 2 * b;
         double blast = 1;

         double acurr = a;
         double atwox = 2 * a;
         double alast = 1;

         double newval = atwox * acurr - alast;
         alast = acurr;
         acurr = newval;

         newval = btwox * bcurr - blast;
         blast = bcurr;
         bcurr = newval;
         double nminus1 = 1.;
         for (unsigned int i = 1; iend != i; ++i) {
            // integrate using recursion relation
            const double c = coeffs[i];
            const double term2 = (blast - alast) / nminus1;

            newval = atwox * acurr - alast;
            alast = acurr;
            acurr = newval;

            newval = btwox * bcurr - blast;
            blast = bcurr;
            bcurr = newval;

            ++nminus1;
            const double term1 = (bcurr - acurr) / (nminus1 + 1.);
            const double intTn = 0.5 * (term1 - term2);
            sum = fast_fma(intTn, c, sum);
         }
      }
   }

   // take care to multiply with the right factor to account for the mapping to
   // normalised range [-1, 1]
   return halfrange * sum;
}

// Clad does not like std::max and std::min so redefined here for simplicity.
inline double max(double x, double y)
{
   return x >= y ? x : y;
}

inline double min(double x, double y)
{
   return x <= y ? x : y;
}

// The last param should be of type bool but it is not as that causes some issues with Cling for some reason...
inline double
poissonIntegral(int code, double mu, double x, double integrandMin, double integrandMax, unsigned int protectNegative)
{
   if (protectNegative && mu < 0.0) {
      return std::exp(-2.0 * mu); // make it fall quickly
   }

   if (code == 1) {
      // Implement integral over x as summation. Add special handling in case
      // range boundaries are not on integer values of x
      integrandMin = max(0, integrandMin);

      if (integrandMax < 0. || integrandMax < integrandMin) {
         return 0;
      }
      const double delta = 100.0 * std::sqrt(mu);
      // If the limits are more than many standard deviations away from the mean,
      // we might as well return the integral of the full Poisson distribution to
      // save computing time.
      if (integrandMin < max(mu - delta, 0.0) && integrandMax > mu + delta) {
         return 1.;
      }

      // The range as integers. ixMin is included, ixMax outside.
      const unsigned int ixMin = integrandMin;
      const unsigned int ixMax = min(integrandMax + 1, (double)std::numeric_limits<unsigned int>::max());

      // Sum from 0 to just before the bin outside of the range.
      if (ixMin == 0) {
         return ROOT::Math::gamma_cdf_c(mu, ixMax, 1);
      } else {
         // If necessary, subtract from 0 to the beginning of the range
         if (ixMin <= mu) {
            return ROOT::Math::gamma_cdf_c(mu, ixMax, 1) - ROOT::Math::gamma_cdf_c(mu, ixMin, 1);
         } else {
            // Avoid catastrophic cancellation in the high tails:
            return ROOT::Math::gamma_cdf(mu, ixMin, 1) - ROOT::Math::gamma_cdf(mu, ixMax, 1);
         }
      }
   }

   // the integral with respect to the mean is the integral of a gamma distribution
   // negative ix does not need protection (gamma returns 0.0)
   const double ix = 1 + x;

   return ROOT::Math::gamma_cdf(integrandMax, ix, 1.0) - ROOT::Math::gamma_cdf(integrandMin, ix, 1.0);
}

inline double logNormalIntegral(double xMin, double xMax, double m0, double k)
{
   const double root2 = std::sqrt(2.);

   double ln_k = TMath::Abs(std::log(k));
   double ret =
      0.5 * (TMath::Erf(std::log(xMax / m0) / (root2 * ln_k)) - TMath::Erf(std::log(xMin / m0) / (root2 * ln_k)));

   return ret;
}

inline double logNormalIntegralStandard(double xMin, double xMax, double mu, double sigma)
{
   const double root2 = std::sqrt(2.);

   double ln_k = std::abs(sigma);
   double ret =
      0.5 * (TMath::Erf((std::log(xMax) - mu) / (root2 * ln_k)) - TMath::Erf((std::log(xMin) - mu) / (root2 * ln_k)));

   return ret;
}

} // namespace RooFitMath

#endif
