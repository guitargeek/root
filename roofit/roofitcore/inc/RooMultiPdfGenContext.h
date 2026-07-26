/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * Authors:                                                                  *
 *   Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MULTIPDF_GEN_CONTEXT
#define ROO_MULTIPDF_GEN_CONTEXT

#include "RooAbsGenContext.h"
#include "RooArgSet.h"

#include <memory>
#include <vector>

class RooMultiPdf;
class RooDataSet;
class RooAbsCategoryLValue;

class RooMultiPdfGenContext : public RooAbsGenContext {
public:
   RooMultiPdfGenContext(const RooMultiPdf &model, const RooArgSet &vars, const RooDataSet *prototype = nullptr,
                         const RooArgSet *auxProto = nullptr, bool verbose = false);

   void setProtoDataOrder(Int_t *lut) override;
   void attach(const RooArgSet &params) override;
   void printMultiline(std::ostream &os, Int_t content, bool verbose = false, TString indent = "") const override;

protected:
   void initGenerator(const RooArgSet &theEvent) override;
   void generateEvent(RooArgSet &theEvent, Int_t remaining) override;

   RooDataSet *createDataSet(const char *name, const char *title, const RooArgSet &obs) override;
   void updateFractions();

   RooAbsCategoryLValue *_idxCat = nullptr;   ///< Clone of index category
   std::unique_ptr<RooArgSet> _idxCatSet;     ///< Owner of index category components
   const RooMultiPdf *_pdf = nullptr;         ///< Original PDF
   std::vector<std::unique_ptr<RooAbsGenContext>> _gcList; ///< List of component generator contexts
   std::vector<int> _gcIndex;                 ///< Index value corresponding to each component
   bool _haveIdxProto = false;                ///< Flag set if generation of index is requested
   TString _idxCatName;                       ///< Name of index category
   Int_t _numPdf = 0;                         ///< Number of generated PDFs
   std::vector<double> _fracThresh;           ///< Fraction threshold array
   std::unique_ptr<RooDataSet> _protoData;    ///<! Prototype dataset

   RooArgSet _allVarsPdf; ///< All pdf variables

   ClassDefOverride(RooMultiPdfGenContext, 0) // Context for generating a dataset from a RooMultiPdf
};

#endif
