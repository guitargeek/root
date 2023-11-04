/*
 * Project: RooFit
 *
 * Copyright (c) 2023, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_RooTable_h
#define RooFit_RooTable_h

#include <RooPrintable.h>

#include <TNamed.h>

/**
\file RooTable.h
\class RooTable
\ingroup Roofitcore

Abstract interface for table objects.
Table objects are the category equivalent of RooPlot objects
(which are used for real-valued objects).
**/

class RooAbsCategory;

/**
\class RooTable
\ingroup Roofitcore

Abstract interface for table objects.
Table objects are the category equivalent of RooPlot objects
(which are used for real-valued objects)
**/

class RooTable : public TNamed, public RooPrintable {
public:
   RooTable() {}
   RooTable(const char *name, const char *title) : TNamed{name, title} {}
   RooTable(const RooTable &other) = default;

   virtual void fill(RooAbsCategory &cat, double weight = 1.0) = 0;

   virtual bool isIdentical(const RooTable &other, bool verbose) = 0;

protected:
   ClassDefOverride(RooTable, 1) // Abstract interface for tables
};

#endif
