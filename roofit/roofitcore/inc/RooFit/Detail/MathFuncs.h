/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2024
 *   Garima Singh, CERN 2023
 *
 * Copyright (c) 2024, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_Detail_MathFuncs_h
#define RooFit_Detail_MathFuncs_h

#include <TMath.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <cmath>

namespace RooFit {

namespace Detail {

namespace MathFuncs {

/// Calculates the binomial coefficient n over k.
/// Equivalent to TMath::Binomial, but inlined.
template <typename RealVal_t>
RealVal_t binomial(int n, int k)
{
   if (n < 0 || k < 0 || n < k)
      return TMath::SignalingNaN();
   if (k == 0 || n == k)
      return 1;

   int k1 = std::min(k, n - k);
   int k2 = n - k1;
   RealVal_t fact = k2 + 1;
   for (RealVal_t i = k1; i > 1.; --i) {
      fact *= (k2 + i) / i;
   }
   return fact;
}

/// The caller needs to make sure that there is at least one coefficient.
template <typename RealVal_t>
RealVal_t bernstein(RealVal_t x, RealVal_t xmin, RealVal_t xmax, RealVal_t *coefs, int nCoefs)
{
   RealVal_t xScaled = (x - xmin) / (xmax - xmin); // rescale to [0,1]
   int degree = nCoefs - 1;                        // n+1 polys of degree n

   // in case list of arguments passed is empty
   if (degree < 0) {
      return TMath::SignalingNaN();
   } else if (degree == 0) {
      return coefs[0];
   } else if (degree == 1) {

      RealVal_t a0 = coefs[0];      // c0
      RealVal_t a1 = coefs[1] - a0; // c1 - c0
      return a1 * xScaled + a0;

   } else if (degree == 2) {

      RealVal_t a0 = coefs[0];            // c0
      RealVal_t a1 = 2 * (coefs[1] - a0); // 2 * (c1 - c0)
      RealVal_t a2 = coefs[2] - a1 - a0;  // c0 - 2 * c1 + c2
      return (a2 * xScaled + a1) * xScaled + a0;
   }

   RealVal_t t = xScaled;
   RealVal_t s = 1. - xScaled;

   RealVal_t result = coefs[0] * s;
   for (int i = 1; i < degree; i++) {
      result = (result + t * binomial<RealVal_t>(degree, i) * coefs[i]) * s;
      t *= xScaled;
   }
   result += t * coefs[degree];

   return result;
}

/// @brief Function to evaluate an un-normalized RooGaussian.
template <typename RealVal_t>
RealVal_t gaussian(RealVal_t x, RealVal_t mean, RealVal_t sigma)
{
   const RealVal_t arg = x - mean;
   const RealVal_t sig = sigma;
   return std::exp(-0.5 * arg * arg / (sig * sig));
}

// RooRatio evaluate function.
template <typename RealVal_t>
RealVal_t ratio(RealVal_t numerator, RealVal_t denominator)
{
   return numerator / denominator;
}

template <typename RealVal_t>
RealVal_t bifurGauss(RealVal_t x, RealVal_t mean, RealVal_t sigmaL, RealVal_t sigmaR)
{
   // Note: this simplification does not work with Clad as of v1.1!
   // return gaussian(x, mean, x < mean ? sigmaL : sigmaR);
   if (x < mean)
      return gaussian(x, mean, sigmaL);
   return gaussian(x, mean, sigmaR);
}

template <typename RealVal_t>
RealVal_t efficiency(RealVal_t effFuncVal, int catIndex, int sigCatIndex)
{
   // Truncate efficiency function in range 0.0-1.0
   effFuncVal = std::clamp(effFuncVal, 0.0, 1.0);

   if (catIndex == sigCatIndex)
      return effFuncVal; // Accept case
   else
      return 1 - effFuncVal; // Reject case
}

/// In pdfMode, a coefficient for the constant term of 1.0 is implied if lowestOrder > 0.
template <bool pdfMode, typename RealVal_t>
RealVal_t polynomial(RealVal_t const *coeffs, int nCoeffs, int lowestOrder, RealVal_t x)
{
   RealVal_t retVal = coeffs[nCoeffs - 1];
   for (int i = nCoeffs - 2; i >= 0; i--)
      retVal = coeffs[i] + x * retVal;
   retVal = retVal * std::pow(x, lowestOrder);
   return retVal + (pdfMode && lowestOrder > 0 ? 1.0 : 0.0);
}

template <typename RealVal_t>
RealVal_t chebychev(RealVal_t *coeffs, unsigned int nCoeffs, RealVal_t x_in, RealVal_t xMin, RealVal_t xMax)
{
   // transform to range [-1, +1]
   const RealVal_t xPrime = (x_in - 0.5 * (xMax + xMin)) / (0.5 * (xMax - xMin));

   // extract current values of coefficients
   RealVal_t sum = 1.;
   if (nCoeffs > 0) {
      RealVal_t curr = xPrime;
      RealVal_t twox = 2 * xPrime;
      RealVal_t last = 1;
      RealVal_t newval = twox * curr - last;
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

template <typename RealVal_t>
RealVal_t constraintSum(RealVal_t const *comp, unsigned int compSize)
{
   RealVal_t sum = 0;
   for (unsigned int i = 0; i < compSize; i++) {
      sum -= std::log(comp[i]);
   }
   return sum;
}

template <typename RealVal_t>
unsigned int getUniformBinning(RealVal_t low, RealVal_t high, RealVal_t val, unsigned int numBins)
{
   RealVal_t binWidth = (high - low) / numBins;
   return val >= high ? numBins - 1 : std::abs((val - low) / binWidth);
}

template <typename RealVal_t>
RealVal_t poisson(RealVal_t x, RealVal_t par)
{
   if (par < 0)
      return TMath::QuietNaN();

   if (x < 0) {
      return 0;
   } else if (x == 0.0) {
      return std::exp(-par);
   } else {
      RealVal_t out = x * std::log(par) - TMath::LnGamma(x + 1.) - par;
      return std::exp(out);
   }
}

template <typename RealVal_t>
RealVal_t flexibleInterpSingle(unsigned int code, RealVal_t low, RealVal_t high, RealVal_t boundary, RealVal_t nominal,
                               RealVal_t paramVal, RealVal_t res)
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
      RealVal_t a = 0.5 * (high + low) - nominal;
      RealVal_t b = 0.5 * (high - low);
      RealVal_t c = 0;
      if (paramVal > 1) {
         return (2 * a + b) * (paramVal - 1) + high - nominal;
      } else if (paramVal < -1) {
         return -1 * (2 * a - b) * (paramVal + 1) + low - nominal;
      } else {
         return a * std::pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 3) {
      // parabolic version of log-normal
      RealVal_t a = 0.5 * (high + low) - nominal;
      RealVal_t b = 0.5 * (high - low);
      RealVal_t c = 0;
      if (paramVal > 1) {
         return (2 * a + b) * (paramVal - 1) + high - nominal;
      } else if (paramVal < -1) {
         return -1 * (2 * a - b) * (paramVal + 1) + low - nominal;
      } else {
         return a * std::pow(paramVal, 2) + b * paramVal + c;
      }
   } else if (code == 4) {
      RealVal_t x = paramVal;
      if (x >= boundary) {
         return x * (high - nominal);
      } else if (x <= -boundary) {
         return x * (nominal - low);
      }

      // interpolate 6th degree
      RealVal_t t = x / boundary;
      RealVal_t eps_plus = high - nominal;
      RealVal_t eps_minus = nominal - low;
      RealVal_t S = 0.5 * (eps_plus + eps_minus);
      RealVal_t A = 0.0625 * (eps_plus - eps_minus);

      return x * (S + t * A * (15 + t * t * (-10 + t * t * 3)));
   } else if (code == 5) {
      RealVal_t x = paramVal;
      RealVal_t mod = 1.0;
      if (x >= boundary) {
         mod = std::pow(high / nominal, +paramVal);
      } else if (x <= -boundary) {
         mod = std::pow(low / nominal, -paramVal);
      } else {
         // interpolate 6th degree exp
         RealVal_t x0 = boundary;

         // GHL: Swagato's suggestions
         RealVal_t powUp = std::pow(high / nominal, x0);
         RealVal_t powDown = std::pow(low / nominal, x0);
         RealVal_t logHi = std::log(high);
         RealVal_t logLo = std::log(low);
         RealVal_t powUpLog = high <= 0.0 ? 0.0 : powUp * logHi;
         RealVal_t powDownLog = low <= 0.0 ? 0.0 : -powDown * logLo;
         RealVal_t powUpLog2 = high <= 0.0 ? 0.0 : powUpLog * logHi;
         RealVal_t powDownLog2 = low <= 0.0 ? 0.0 : -powDownLog * logLo;

         RealVal_t S0 = 0.5 * (powUp + powDown);
         RealVal_t A0 = 0.5 * (powUp - powDown);
         RealVal_t S1 = 0.5 * (powUpLog + powDownLog);
         RealVal_t A1 = 0.5 * (powUpLog - powDownLog);
         RealVal_t S2 = 0.5 * (powUpLog2 + powDownLog2);
         RealVal_t A2 = 0.5 * (powUpLog2 - powDownLog2);

         // fcns+der+2nd_der are eq at bd

         RealVal_t a = 1. / (8 * x0) * (15 * A0 - 7 * x0 * S1 + x0 * x0 * A2);
         RealVal_t b = 1. / (8 * x0 * x0) * (-24 + 24 * S0 - 9 * x0 * A1 + x0 * x0 * S2);
         RealVal_t c = 1. / (4 * std::pow(x0, 3)) * (-5 * A0 + 5 * x0 * S1 - x0 * x0 * A2);
         RealVal_t d = 1. / (4 * std::pow(x0, 4)) * (12 - 12 * S0 + 7 * x0 * A1 - x0 * x0 * S2);
         RealVal_t e = 1. / (8 * std::pow(x0, 5)) * (+3 * A0 - 3 * x0 * S1 + x0 * x0 * A2);
         RealVal_t f = 1. / (8 * std::pow(x0, 6)) * (-8 + 8 * S0 - 5 * x0 * A1 + x0 * x0 * S2);

         // evaluate the 6-th degree polynomial using Horner's method
         RealVal_t value = 1. + x * (a + x * (b + x * (c + x * (d + x * (e + x * f)))));
         mod = value;
      }
      return res * (mod - 1.0);
   }

   return 0.0;
}

template <typename RealVal_t>
RealVal_t flexibleInterp(unsigned int code, RealVal_t const *params, unsigned int n, RealVal_t const *low,
                         RealVal_t const *high, RealVal_t boundary, RealVal_t nominal, int doCutoff)
{
   RealVal_t total = nominal;
   for (std::size_t i = 0; i < n; ++i) {
      total += flexibleInterpSingle(code, low[i], high[i], boundary, nominal, params[i], total);
   }

   return doCutoff && total <= 0 ? TMath::Limits<RealVal_t>::Min() : total;
}

template <typename RealVal_t>
RealVal_t landau(RealVal_t x, RealVal_t mu, RealVal_t sigma)
{
   if (sigma <= 0.)
      return 0.;
   return ROOT::Math::landau_pdf((x - mu) / sigma);
}

template <typename RealVal_t>
RealVal_t logNormal(RealVal_t x, RealVal_t k, RealVal_t m0)
{
   return ROOT::Math::lognormal_pdf(x, std::log(m0), std::abs(std::log(k)));
}

template <typename RealVal_t>
RealVal_t logNormalStandard(RealVal_t x, RealVal_t sigma, RealVal_t mu)
{
   return ROOT::Math::lognormal_pdf(x, mu, std::abs(sigma));
}

template <typename RealVal_t>
RealVal_t effProd(RealVal_t eff, RealVal_t pdf)
{
   return eff * pdf;
}

template <typename RealVal_t>
RealVal_t nll(RealVal_t pdf, RealVal_t weight, int binnedL, int doBinOffset)
{
   if (binnedL) {
      // Special handling of this case since std::log(Poisson(0,0)=0 but can't be
      // calculated with usual log-formula since std::log(mu)=0. No update of result
      // is required since term=0.
      if (std::abs(pdf) < 1e-10 && std::abs(weight) < 1e-10) {
         return 0.0;
      }
      if (doBinOffset) {
         return pdf - weight - weight * (std::log(pdf) - std::log(weight));
      }
      return pdf - weight * std::log(pdf) + TMath::LnGamma(weight + 1);
   } else {
      return -weight * std::log(pdf);
   }
}

template <typename RealVal_t>
RealVal_t recursiveFraction(RealVal_t *a, unsigned int n)
{
   RealVal_t prod = a[0];

   for (unsigned int i = 1; i < n; ++i) {
      prod *= 1.0 - a[i];
   }

   return prod;
}

template <typename RealVal_t>
RealVal_t cbShape(RealVal_t m, RealVal_t m0, RealVal_t sigma, RealVal_t alpha, RealVal_t n)
{
   RealVal_t t = (m - m0) / sigma;
   if (alpha < 0)
      t = -t;

   RealVal_t absAlpha = std::abs((RealVal_t)alpha);

   if (t >= -absAlpha) {
      return std::exp(-0.5 * t * t);
   } else {
      RealVal_t a = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
      RealVal_t b = n / absAlpha - absAlpha;

      return a / std::pow(b - t, n);
   }
}

// For RooCBShape
template <typename RealVal_t>
RealVal_t approxErf(RealVal_t arg)
{
   if (arg > 5.0)
      return 1.0;
   if (arg < -5.0)
      return -1.0;

   return std::erf(arg);
}

/// @brief Function to calculate the integral of an un-normalized RooGaussian over x. To calculate the integral over
/// mean, just interchange the respective values of x and mean.
/// @param xMin Minimum value of variable to integrate wrt.
/// @param xMax Maximum value of of variable to integrate wrt.
/// @param mean Mean.
/// @param sigma Sigma.
/// @return The integral of an un-normalized RooGaussian over the value in x.
template <typename RealVal_t>
RealVal_t gaussianIntegral(RealVal_t xMin, RealVal_t xMax, RealVal_t mean, RealVal_t sigma)
{
   // The normalisation constant 1./sqrt(2*pi*sigma^2) is left out in evaluate().
   // Therefore, the integral is scaled up by that amount to make RooFit normalise
   // correctly.
   RealVal_t resultScale = 0.5 * std::sqrt(TMath::TwoPi()) * sigma;

   // Here everything is scaled and shifted into a standard normal distribution:
   RealVal_t xscale = TMath::Sqrt2() * sigma;
   RealVal_t scaledMin = 0.;
   RealVal_t scaledMax = 0.;
   scaledMin = (xMin - mean) / xscale;
   scaledMax = (xMax - mean) / xscale;

   // Here we go for maximum precision: We compute all integrals in the UPPER
   // tail of the Gaussian, because erfc has the highest precision there.
   // Therefore, the different cases for range limits in the negative hemisphere are mapped onto
   // the equivalent points in the upper hemisphere using erfc(-x) = 2. - erfc(x)
   RealVal_t ecmin = std::erfc(std::abs(scaledMin));
   RealVal_t ecmax = std::erfc(std::abs(scaledMax));

   RealVal_t cond = 0.0;
   // Don't put this "prd" inside the "if" because clad will not be able to differentiate the code correctly (as of
   // v1.1)!
   RealVal_t prd = scaledMin * scaledMax;
   if (prd < 0.0) {
      cond = 2.0 - (ecmin + ecmax);
   } else if (scaledMax <= 0.0) {
      cond = ecmax - ecmin;
   } else {
      cond = ecmin - ecmax;
   }
   return resultScale * cond;
}

template <typename RealVal_t>
RealVal_t bifurGaussIntegral(RealVal_t xMin, RealVal_t xMax, RealVal_t mean, RealVal_t sigmaL, RealVal_t sigmaR)
{
   const RealVal_t xscaleL = TMath::Sqrt2() * sigmaL;
   const RealVal_t xscaleR = TMath::Sqrt2() * sigmaR;

   const RealVal_t resultScale = 0.5 * std::sqrt(TMath::TwoPi());

   if (xMax < mean) {
      return resultScale * (sigmaL * (std::erf((xMax - mean) / xscaleL) - std::erf((xMin - mean) / xscaleL)));
   } else if (xMin > mean) {
      return resultScale * (sigmaR * (std::erf((xMax - mean) / xscaleR) - std::erf((xMin - mean) / xscaleR)));
   } else {
      return resultScale *
             (sigmaR * std::erf((xMax - mean) / xscaleR) - sigmaL * std::erf((xMin - mean) / xscaleL));
   }
}

template <typename RealVal_t>
RealVal_t exponentialIntegral(RealVal_t xMin, RealVal_t xMax, RealVal_t constant)
{
   if (constant == 0.0) {
      return xMax - xMin;
   }

   return (std::exp(constant * xMax) - std::exp(constant * xMin)) / constant;
}

/// In pdfMode, a coefficient for the constant term of 1.0 is implied if lowestOrder > 0.
template <bool pdfMode, typename RealVal_t>
RealVal_t polynomialIntegral(RealVal_t const *coeffs, int nCoeffs, int lowestOrder, RealVal_t xMin, RealVal_t xMax)
{
   int denom = lowestOrder + nCoeffs;
   RealVal_t min = coeffs[nCoeffs - 1] / RealVal_t(denom);
   RealVal_t max = coeffs[nCoeffs - 1] / RealVal_t(denom);

   for (int i = nCoeffs - 2; i >= 0; i--) {
      denom--;
      min = (coeffs[i] / RealVal_t(denom)) + xMin * min;
      max = (coeffs[i] / RealVal_t(denom)) + xMax * max;
   }

   max = max * std::pow(xMax, 1 + lowestOrder);
   min = min * std::pow(xMin, 1 + lowestOrder);

   return max - min + (pdfMode && lowestOrder > 0.0 ? xMax - xMin : 0.0);
}

/// use fast FMA if available, fall back to normal arithmetic if not
template <typename RealVal_t>
RealVal_t fast_fma(RealVal_t x, RealVal_t y, RealVal_t z) noexcept
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

template <typename RealVal_t>
RealVal_t chebychevIntegral(RealVal_t const *coeffs, unsigned int nCoeffs, RealVal_t xMin, RealVal_t xMax,
                            RealVal_t xMinFull, RealVal_t xMaxFull)
{
   const RealVal_t halfrange = .5 * (xMax - xMin);
   const RealVal_t mid = .5 * (xMax + xMin);

   // the full range of the function is mapped to the normalised [-1, 1] range
   const RealVal_t b = (xMaxFull - mid) / halfrange;
   const RealVal_t a = (xMinFull - mid) / halfrange;

   // coefficient for integral(T_0(x)) is 1 (implicit), integrate by hand
   // T_0(x) and T_1(x), and use for n > 1: integral(T_n(x) dx) =
   // (T_n+1(x) / (n + 1) - T_n-1(x) / (n - 1)) / 2
   RealVal_t sum = b - a; // integrate T_0(x) by hand

   const unsigned int iend = nCoeffs;
   if (iend > 0) {
      {
         // integrate T_1(x) by hand...
         const RealVal_t c = coeffs[0];
         sum = fast_fma(0.5 * (b + a) * (b - a), c, sum);
      }
      if (1 < iend) {
         RealVal_t bcurr = b;
         RealVal_t btwox = 2 * b;
         RealVal_t blast = 1;

         RealVal_t acurr = a;
         RealVal_t atwox = 2 * a;
         RealVal_t alast = 1;

         RealVal_t newval = atwox * acurr - alast;
         alast = acurr;
         acurr = newval;

         newval = btwox * bcurr - blast;
         blast = bcurr;
         bcurr = newval;
         RealVal_t nminus1 = 1.;
         for (unsigned int i = 1; iend != i; ++i) {
            // integrate using recursion relation
            const RealVal_t c = coeffs[i];
            const RealVal_t term2 = (blast - alast) / nminus1;

            newval = atwox * acurr - alast;
            alast = acurr;
            acurr = newval;

            newval = btwox * bcurr - blast;
            blast = bcurr;
            bcurr = newval;

            ++nminus1;
            const RealVal_t term1 = (bcurr - acurr) / (nminus1 + 1.);
            const RealVal_t intTn = 0.5 * (term1 - term2);
            sum = fast_fma(intTn, c, sum);
         }
      }
   }

   // take care to multiply with the right factor to account for the mapping to
   // normalised range [-1, 1]
   return halfrange * sum;
}

// Clad does not like std::max and std::min so redefined here for simplicity.
template <typename RealVal_t>
RealVal_t max(RealVal_t x, RealVal_t y)
{
   return x >= y ? x : y;
}

template <typename RealVal_t>
RealVal_t min(RealVal_t x, RealVal_t y)
{
   return x <= y ? x : y;
}

// The last param should be of type bool but it is not as that causes some issues with Cling for some reason...
template <typename RealVal_t>
RealVal_t poissonIntegral(int code, RealVal_t mu, RealVal_t x, RealVal_t integrandMin, RealVal_t integrandMax,
                          unsigned int protectNegative)
{
   if (protectNegative && mu < 0.0) {
      return std::exp(-2.0 * mu); // make it fall quickly
   }

   if (code == 1) {
      // Implement integral over x as summation. Add special handling in case
      // range boundaries are not on integer values of x
      integrandMin = max<RealVal_t>(0., integrandMin);

      if (integrandMax < 0. || integrandMax < integrandMin) {
         return 0;
      }
      const RealVal_t delta = 100.0 * std::sqrt(mu);
      // If the limits are more than many standard deviations away from the mean,
      // we might as well return the integral of the full Poisson distribution to
      // save computing time.
      if (integrandMin < max<RealVal_t>(mu - delta, 0.0) && integrandMax > mu + delta) {
         return 1.;
      }

      // The range as integers. ixMin is included, ixMax outside.
      const unsigned int ixMin = integrandMin;
      const unsigned int ixMax = min<RealVal_t>(integrandMax + 1, (RealVal_t)std::numeric_limits<unsigned int>::max());

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
   const RealVal_t ix = 1 + x;

   return ROOT::Math::gamma_cdf(integrandMax, ix, 1.0) - ROOT::Math::gamma_cdf(integrandMin, ix, 1.0);
}

template <typename RealVal_t>
RealVal_t logNormalIntegral(RealVal_t xMin, RealVal_t xMax, RealVal_t m0, RealVal_t k)
{
   const RealVal_t root2 = std::sqrt(2.);

   RealVal_t ln_k = std::abs(std::log(k));
   RealVal_t ret =
      0.5 * (std::erf(std::log(xMax / m0) / (root2 * ln_k)) - std::erf(std::log(xMin / m0) / (root2 * ln_k)));

   return ret;
}

template <typename RealVal_t>
RealVal_t logNormalIntegralStandard(RealVal_t xMin, RealVal_t xMax, RealVal_t mu, RealVal_t sigma)
{
   const RealVal_t root2 = std::sqrt(2.);

   RealVal_t ln_k = std::abs(sigma);
   RealVal_t ret =
      0.5 * (std::erf((std::log(xMax) - mu) / (root2 * ln_k)) - std::erf((std::log(xMin) - mu) / (root2 * ln_k)));

   return ret;
}

template <typename RealVal_t>
RealVal_t cbShapeIntegral(RealVal_t mMin, RealVal_t mMax, RealVal_t m0, RealVal_t sigma, RealVal_t alpha, RealVal_t n)
{
   const RealVal_t sqrtPiOver2 = 1.2533141373;
   const RealVal_t sqrt2 = 1.4142135624;

   RealVal_t result = 0.0;
   bool useLog = false;

   if (std::abs(n - 1.0) < 1.0e-05)
      useLog = true;

   RealVal_t sig = std::abs(sigma);

   RealVal_t tmin = (mMin - m0) / sig;
   RealVal_t tmax = (mMax - m0) / sig;

   if (alpha < 0) {
      RealVal_t tmp = tmin;
      tmin = -tmax;
      tmax = -tmp;
   }

   RealVal_t absAlpha = std::abs(alpha);

   if (tmin >= -absAlpha) {
      result += sig * sqrtPiOver2 * (approxErf(tmax / sqrt2) - approxErf(tmin / sqrt2));
   } else if (tmax <= -absAlpha) {
      RealVal_t a = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
      RealVal_t b = n / absAlpha - absAlpha;

      if (useLog) {
         result += a * sig * (std::log(b - tmin) - std::log(b - tmax));
      } else {
         result += a * sig / (1.0 - n) * (1.0 / (std::pow(b - tmin, n - 1.0)) - 1.0 / (std::pow(b - tmax, n - 1.0)));
      }
   } else {
      RealVal_t a = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
      RealVal_t b = n / absAlpha - absAlpha;

      RealVal_t term1 = 0.0;
      if (useLog) {
         term1 = a * sig * (std::log(b - tmin) - std::log(n / absAlpha));
      } else {
         term1 = a * sig / (1.0 - n) * (1.0 / (std::pow(b - tmin, n - 1.0)) - 1.0 / (std::pow(n / absAlpha, n - 1.0)));
      }

      RealVal_t term2 = sig * sqrtPiOver2 * (approxErf(tmax / sqrt2) - approxErf(-absAlpha / sqrt2));

      result += term1 + term2;
   }

   if (result == 0)
      return 1.E-300;
   return result;
}

template <typename RealVal_t>
RealVal_t bernsteinIntegral(RealVal_t xlo, RealVal_t xhi, RealVal_t xmin, RealVal_t xmax, RealVal_t *coefs, int nCoefs)
{
   RealVal_t xloScaled = (xlo - xmin) / (xmax - xmin);
   RealVal_t xhiScaled = (xhi - xmin) / (xmax - xmin);

   int degree = nCoefs - 1; // n+1 polys of degree n
   RealVal_t norm = 0.;

   for (int i = 0; i <= degree; ++i) {
      // for each of the i Bernstein basis polynomials
      // represent it in the 'power basis' (the naive polynomial basis)
      // where the integral is straight forward.
      RealVal_t temp = 0.;
      for (int j = i; j <= degree; ++j) { // power basis≈ß
         RealVal_t binCoefs = binomial<RealVal_t>(degree, j) * binomial<RealVal_t>(j, i);
         RealVal_t oneOverJPlusOne = 1. / (j + 1.);
         RealVal_t powDiff = std::pow(xhiScaled, j + 1.) - std::pow(xloScaled, j + 1.);
         temp += std::pow(-1., j - i) * binCoefs * powDiff * oneOverJPlusOne;
      }
      temp *= coefs[i]; // include coeff
      norm += temp;     // add this basis's contribution to total
   }

   return norm * (xmax - xmin);
}

} // namespace MathFuncs

} // namespace Detail

} // namespace RooFit

#endif
