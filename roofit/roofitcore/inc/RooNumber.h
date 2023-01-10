/*
 * Project: RooFit
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file RooNumber.h
\class RooNumber
\ingroup Roofitcore

RooNumber implements numeric constants used by RooFit.
**/

#ifndef RooFit_RooNumber_h
#define RooFit_RooNumber_h

class RooNumber {
public:
   /// Return internal infinity representation.
   static constexpr double infinity() { return _inf; }

   /// Return sign(x) if x is infinite by the RooFit-internal specification,
   /// otherwise return zero.
   static constexpr int isInfinite(double x) { return (x >= +_inf) ? +1 : ((x <= -_inf) ? -1 : 0); }

private:
   // This assumes a well behaved IEEE-754 floating point implementation. The
   // next line may generate a compiler warning that can be ignored.
   static constexpr double _inf = 1.e30;
};

#endif
