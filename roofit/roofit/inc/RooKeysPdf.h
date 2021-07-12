/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooKeysPdf.h,v 1.10 2007/05/11 09:13:07 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   UC San Diego,        raven@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_KEYS
#define ROO_KEYS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooKeysPdf : public RooAbsPdf {
public:
  enum Mirror { NoMirror, MirrorLeft, MirrorRight, MirrorBoth,
      MirrorAsymLeft, MirrorAsymLeftRight,
      MirrorAsymRight, MirrorLeftAsymRight,
      MirrorAsymBoth };
  RooKeysPdf() ;
  RooKeysPdf(const char *name, const char *title,
             RooAbsReal& x, RooDataSet& data, Mirror mirror= NoMirror,
        Double_t rho=1);
  RooKeysPdf(const char *name, const char *title,
             RooAbsReal& x, RooRealVar& xdata, RooDataSet& data, Mirror mirror= NoMirror,
        Double_t rho=1);
  RooKeysPdf(const char *name, const char *title, RooAbsRealLValue& x, double const* lookupTable);
  RooKeysPdf(const RooKeysPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const {return new RooKeysPdf(*this,newname); }
  virtual ~RooKeysPdf();

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
     const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
  virtual Int_t getMaxVal(const RooArgSet& vars) const;
  virtual Double_t maxVal(Int_t code) const;

  struct Configuration {
    bool mirrorLeft = false;
    bool mirrorRight = false;
    bool asymLeft = false;
    bool asymRight = false;
  };

  void LoadDataSet( RooDataSet& data, const char* xName, double rho, RooKeysPdf::Configuration const& cfg);

  /// Returns pointer to the beginning of the lookup table that defines this RooKeysPdf.
  double const* lookupTable() const { return _lookupTable; }
  /// Returns the number of points used for the lookup table of this RooKeysPdf.
  int nPoints() const { return _nPoints; }

protected:

  RooRealProxy _x ;
  Double_t evaluate() const;

private:
  // how far you have to go out in a Gaussian until it is smaller than the
  // machine precision
  static const Double_t _nSigma; //!

  enum { _nPoints = 1000 };
  double _lookupTable[_nPoints+1];

  double g(int nEvents, std::vector<double> const& dataPts, double x, double sigma) const;

  // cached info on variable
  double _lo;
  double _hi;
  double _binWidth;

  ClassDef(RooKeysPdf,3) // One-dimensional non-parametric kernel estimation p.d.f.
};

#endif
