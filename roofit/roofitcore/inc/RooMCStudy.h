/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooMCStudy.h,v 1.18 2007/05/11 10:14:56 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MC_STUDY
#define ROO_MC_STUDY

#include "TList.h"
#include "TNamed.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooFit/Kwargs.h"
#include <list>
class RooAbsPdf;
class RooDataSet ;
class RooAbsData ;
class RooAbsGenContext ;
class RooFitResult ;
class RooPlot ;
class RooRealVar ;
class RooAbsMCStudyModule ;

class RooMCStudy : public TNamed {
public:
   /// Keyword arguments for the RooMCStudy constructor.
   /// Each member corresponds to the RooFit command argument of the same name.
   struct Kwargs {
      RooFit::Detail::Kwarg::Binned Binned;
      RooFit::Detail::Kwarg::ConditionalObservables ConditionalObservables;
      RooFit::Detail::Kwarg::Constrain Constrain;
      RooFit::Detail::Kwarg::Extended Extended;
      RooFit::Detail::Kwarg::ExternalConstraints ExternalConstraints;
      RooFit::Detail::Kwarg::FitModel FitModel;
      RooFit::Detail::Kwarg::FitOptions FitOptions;
      RooFit::Detail::Kwarg::ProjectedObservables ProjectedObservables;
      RooFit::Detail::Kwarg::ProtoData ProtoData;
      RooFit::Detail::Kwarg::Silence Silence;
      RooFit::Detail::Kwarg::Verbose Verbose;

      std::vector<RooCmdArg const *> cmdArgs() const
      {
         // clang-format off
         return {&Binned, &ConditionalObservables, &Constrain, &Extended, &ExternalConstraints, &FitModel, &FitOptions,
                 &ProjectedObservables, &ProtoData, &Silence, &Verbose};
         // clang-format on
      }
   };

  RooMCStudy(const RooAbsPdf& model, const RooArgSet& observables,
        const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
             const RooCmdArg& arg3={}, const RooCmdArg& arg4={}, const RooCmdArg& arg5={},
             const RooCmdArg& arg6={}, const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooMCStudy(const RooAbsPdf &model, const RooArgSet &observables, Kwargs const &kwargs, const RooCmdArg &arg1 = {},
             const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
             const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
             const RooCmdArg &arg8 = {})
     : RooMCStudy(RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                                 "RooMCStudy::RooMCStudy()"),
                  model, observables)
  {
  }

  RooMCStudy(const RooMCStudy&) = delete;
  RooMCStudy(RooMCStudy &&) = delete;
  RooMCStudy& operator=(const RooMCStudy&) = delete;
  RooMCStudy& operator=(RooMCStudy &&) = delete;

  ~RooMCStudy() override ;

  // Method to add study modules
  void addModule(RooAbsMCStudyModule& module) ;


  // Run methods
  bool generateAndFit(Int_t nSamples, Int_t nEvtPerSample=0, bool keepGenData=false, const char* asciiFilePat=nullptr) ;
  bool generate(Int_t nSamples, Int_t nEvtPerSample=0, bool keepGenData=false, const char* asciiFilePat=nullptr) ;
  bool fit(Int_t nSamples, const char* asciiFilePat) ;
  bool fit(Int_t nSamples, TList& dataSetList) ;
  bool addFitResult(const RooFitResult& fr) ;

  // Result accessors
  const RooArgSet* fitParams(Int_t sampleNum) const ;
  const RooFitResult* fitResult(Int_t sampleNum) const ;
        RooAbsData* genData(Int_t sampleNum) const ;
  const RooDataSet& fitParDataSet() ;
  /// Return dataset with the generator parameter values used for each toy, including any
  /// modification by constraint-p.d.f. sampling or by study modules. When constraints are
  /// used, the values are sampled from the constraint p.d.f.s for each toy and thus generally
  /// differ from the fitted parameter values. Unlike fitParDataSet(), this dataset has one
  /// entry per generated toy, including toys for which the fit did not converge.
  const RooDataSet* genParDataSet() const {
    return _genParData.get();
  }

  // Plot methods
  /// Keyword arguments for RooMCStudy::plotParam().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct PlotParamKwargs {
     RooFit::Detail::Kwarg::AddTo AddTo;
     RooFit::Detail::Kwarg::Asymmetry Asymmetry;
     RooFit::Detail::Kwarg::Binning Binning;
     RooFit::Detail::Kwarg::Bins Bins;
     RooFit::Detail::Kwarg::Cut Cut;
     RooFit::Detail::Kwarg::CutRange CutRange;
     RooFit::Detail::Kwarg::DataError DataError;
     RooFit::Detail::Kwarg::DrawOption DrawOption;
     RooFit::Detail::Kwarg::Efficiency Efficiency;
     RooFit::Detail::Kwarg::FillColor FillColor;
     RooFit::Detail::Kwarg::FillStyle FillStyle;
     RooFit::Detail::Kwarg::FitGauss FitGauss;
     RooFit::Detail::Kwarg::Frame Frame;
     RooFit::Detail::Kwarg::FrameBins FrameBins;
     RooFit::Detail::Kwarg::FrameRange FrameRange;
     RooFit::Detail::Kwarg::Invisible Invisible;
     RooFit::Detail::Kwarg::LineColor LineColor;
     RooFit::Detail::Kwarg::LineStyle LineStyle;
     RooFit::Detail::Kwarg::LineWidth LineWidth;
     RooFit::Detail::Kwarg::MarkerColor MarkerColor;
     RooFit::Detail::Kwarg::MarkerSize MarkerSize;
     RooFit::Detail::Kwarg::MarkerStyle MarkerStyle;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::Range Range;
     RooFit::Detail::Kwarg::RefreshNorm RefreshNorm;
     RooFit::Detail::Kwarg::Rescale Rescale;
     RooFit::Detail::Kwarg::XErrorSize XErrorSize;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&AddTo, &Asymmetry, &Binning, &Bins, &Cut, &CutRange, &DataError, &DrawOption, &Efficiency, &FillColor,
                &FillStyle, &FitGauss, &Frame, &FrameBins, &FrameRange, &Invisible, &LineColor, &LineStyle, &LineWidth,
                &MarkerColor, &MarkerSize, &MarkerStyle, &Name, &Range, &RefreshNorm, &Rescale, &XErrorSize};
        // clang-format on
     }
  };

  RooPlot* plotParamOn(RooPlot* frame, const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                       const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                       const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                       const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotParamOn(RooPlot *frame, RooAbsData::PlotOnKwargs const &kwargs, const RooCmdArg &arg1 = {},
                       const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                       const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                       const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotParamOn()");
     return plotParamOn(frame, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooPlot* plotParam(const RooRealVar& param, const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                     const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                     const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                     const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotParam(const RooRealVar &param, PlotParamKwargs const &kwargs, const RooCmdArg &arg1 = {},
                     const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                     const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                     const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotParam()");
     return plotParam(param, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooPlot* plotParam(const char* paramName, const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                     const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                     const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                     const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotParam(const char *paramName, PlotParamKwargs const &kwargs, const RooCmdArg &arg1 = {},
                     const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                     const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                     const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotParam()");
     return plotParam(paramName, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooPlot* plotNLL(const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                     const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                     const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                     const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotNLL(PlotParamKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                   const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                   const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotNLL()");
     return plotNLL(*a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooPlot* plotError(const RooRealVar& param, const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                     const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                     const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                     const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotError(const RooRealVar &param, PlotParamKwargs const &kwargs, const RooCmdArg &arg1 = {},
                     const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                     const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                     const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotError()");
     return plotError(param, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooPlot* plotPull(const RooRealVar& param, const RooCmdArg& arg1, const RooCmdArg& arg2={},
                     const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                     const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                     const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *plotPull(const RooRealVar &param, PlotParamKwargs const &kwargs, const RooCmdArg &arg1 = {},
                    const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                    const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                    const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooMCStudy::plotPull()");
     return plotPull(param, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }

  RooPlot* plotNLL(double lo, double hi, Int_t nBins=100) ;
  RooPlot* plotError(const RooRealVar& param, double lo, double hi, Int_t nbins=100) ;
  RooPlot* plotPull(const RooRealVar& param, double lo=-3.0, double hi=3.0, Int_t nbins=25, bool fitGauss=false) ;

protected:
   // Delegated to by the constructor that takes kwargs.
   RooMCStudy(std::array<RooCmdArg const *, 8> const &a, const RooAbsPdf &model, const RooArgSet &observables)
      : RooMCStudy(model, observables, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7])
   {
   }

  friend class RooAbsMCStudyModule ;

  RooPlot* makeFrameAndPlotCmd(const RooRealVar& param, RooLinkedList& cmdList, bool symRange=false) const ;

  bool run(bool generate, bool fit, Int_t nSamples, Int_t nEvtPerSample, bool keepGenData, const char* asciiFilePat) ;
  bool fitSample(RooAbsData* genSample) ;
  RooFit::OwningPtr<RooFitResult> doFit(RooAbsData* genSample) ;

  void calcPulls() ;

  RooAbsData*       _genSample = nullptr;       ///< Currently generated sample
  RooAbsPdf*        _genModel = nullptr;        ///< Generator model
  std::unique_ptr<RooAbsGenContext> _genContext;      ///< Generator context
  RooArgSet        _genInitParams;   ///< List of original generator parameters
  RooArgSet        _genParams;       ///< List of actual generator parameters
  const RooDataSet* _genProtoData ;    ///< Generator prototype data set
  RooArgSet         _projDeps ;        ///< List of projected dependents in fit

  std::unique_ptr<RooAbsPdf>        _constrPdf;        ///< Constraints p.d.f
  std::unique_ptr<RooAbsGenContext> _constrGenContext; ///< Generator context for constraints p.d.f

  RooArgSet    _dependents ;    ///< List of dependents
  RooArgSet    _allDependents ; ///< List of generate + prototype dependents
  RooAbsPdf*   _fitModel ;      ///< Fit model
  RooArgSet   _fitInitParams; ///< List of initial values of fit parameters
  RooArgSet   _fitParams;     ///< List of actual fit parameters
  std::unique_ptr<RooRealVar>  _nllVar;
  std::unique_ptr<RooRealVar>  _ngenVar;

  TList       _genDataList ;    // List of generated data sample
  TList       _fitResList ;     // List of RooFitResult fit output objects
  std::unique_ptr<RooDataSet> _genParData;     // List of generated parameters of each sample
  std::unique_ptr<RooDataSet> _fitParData;     // Data set of fit parameters of each sample
  RooLinkedList _fitOptList ;   // Fit option command list
  bool      _extendedGen ;    // Add poisson term to number of events to generate?
  bool      _binGenData ;     // Bin data between generating and fitting
  double    _nExpGen ;        // Number of expected events to generate in extended mode
  bool      _randProto ;      // Randomize order of prototype data access

  bool      _canAddFitResults ; ///< Allow adding of external fit results?
  bool      _verboseGen       ; ///< Verbose generation?
  bool      _perExptGenParams = false; ///< Do generation parameter change per event?
  bool      _silence          ; ///< Silent running mode?

  std::list<RooAbsMCStudyModule*> _modList ; ///< List of additional study modules ;

  // Utilities for modules ;
  RooFit::OwningPtr<RooFitResult> refit(RooAbsData* genSample=nullptr) ;
  void resetFitParams() ;
  void RecursiveRemove(TObject *obj) override;

  ClassDefOverride(RooMCStudy,0) // A general purpose toy Monte Carlo study manager
} ;


#endif

