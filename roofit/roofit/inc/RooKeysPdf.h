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

#include <vector>

class RooRealVar;

class RooKeysPdf : public RooAbsPdf {
public:
  /// Boundary correction obtained by reflecting the data across the lower
  /// and/or upper edge of the observable range. Symmetric mirroring adds the
  /// reflected events (density flat at the boundary), asymmetric mirroring
  /// subtracts them (density vanishing at the boundary).
  enum Mirror {
     NoMirror,            ///< No boundary correction
     MirrorLeft,          ///< Symmetric mirror at the lower edge
     MirrorRight,         ///< Symmetric mirror at the upper edge
     MirrorBoth,          ///< Symmetric mirror at both edges
     MirrorAsymLeft,      ///< Asymmetric mirror at the lower edge
     MirrorAsymLeftRight, ///< Asymmetric mirror at the lower edge, symmetric mirror at the upper edge
     MirrorAsymRight,     ///< Asymmetric mirror at the upper edge
     MirrorLeftAsymRight, ///< Symmetric mirror at the lower edge, asymmetric mirror at the upper edge
     MirrorAsymBoth       ///< Asymmetric mirror at both edges
  };
  RooKeysPdf() ;
  RooKeysPdf(const char *name, const char *title,
             RooAbsReal& xpdf, RooDataSet& data, Mirror mirror= NoMirror,
        double rho=1);
  RooKeysPdf(const char *name, const char *title,
             RooAbsReal& x, RooRealVar& xdata, RooDataSet& data, Mirror mirror= NoMirror,
        double rho=1);
  RooKeysPdf(const char *name, const char *title, RooAbsRealLValue& x, double const* lookupTable);
  RooKeysPdf(const RooKeysPdf& other, const char* name=nullptr);
  TObject* clone(const char* newname=nullptr) const override {return new RooKeysPdf(*this,newname); }
  ~RooKeysPdf() override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
     const char* rangeName = nullptr) const override;
  double analyticalIntegral(Int_t code, const char* rangeName = nullptr) const override;
  Int_t getMaxVal(const RooArgSet& vars) const override;
  double maxVal(Int_t code) const override;

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
  double evaluate() const override;

private:
  // how far you have to go out in a Gaussian until it is smaller than the
  // machine precision
  static const double _nSigma; ///<!

  constexpr static int _nPoints{1000};
  double _lookupTable[_nPoints+1];

  double g(std::vector<double> const& dataPts, double x, double sigma) const;

  // cached info on variable
  double _lo;
  double _hi;
  double _binWidth;

  ClassDefOverride(RooKeysPdf,3) // One-dimensional non-parametric kernel estimation p.d.f.
};

#endif
