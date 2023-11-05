/*
 * Project: RooFit
 * Authors:
 *   Emmanouil Michalainas, CERN 3 March 2021
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

/**
\file Batches.h
\class Batch
\class Batches
\ingroup roofit_dev_docs_batchcompute

These classes encapsulate the necessary data for the computations.
They are lightweight objects designed to be passed by value and also flexible,
so that they can contain data for every kind of compute function.
**/
#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_H
#define ROOFIT_BATCHCOMPUTE_BATCHES_H

#include <cstdint>

namespace RooBatchCompute {

class Batch {
public:
   const double *__restrict _array = nullptr;
   bool _isVector = false;

   // Access the value for event `i`. Non-vector (scalar) inputs always return
   // their single value, so the same access pattern is safe on CPU and GPU.
   __roodevice__ constexpr double operator[](std::size_t i) const noexcept { return _isVector ? _array[i] : _array[0]; }
};

class Batches {
public:
   Batch *args = nullptr;
   double *extra;
   std::size_t nEvents = 0;
   std::size_t nBatches = 0;
   std::size_t nExtra = 0;
   double *__restrict output = nullptr;
};

} // end namespace RooBatchCompute

#endif // #ifdef ROOFIT_BATCHCOMPUTE_BATCHES_H
