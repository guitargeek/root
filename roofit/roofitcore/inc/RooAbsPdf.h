/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsPdf.h,v 1.90 2007/07/21 21:32:52 wouter Exp $
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
#ifndef ROO_ABS_PDF
#define ROO_ABS_PDF

#include <RooAbsReal.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFit/Kwargs.h>
#include <RooFit/UniqueId.h>
#include <RooGlobalFunc.h>
#include <RooObjCacheManager.h>

class RooArgSet ;
class RooAbsGenContext ;
class RooFitResult ;
class RooNumGenConfig ;


class RooAbsPdf : public RooAbsReal {
public:

  // Constructors, assignment etc
  RooAbsPdf() ;
  RooAbsPdf(const char *name, const char *title=nullptr) ;
  RooAbsPdf(const char *name, const char *title, double minVal, double maxVal) ;
  // RooAbsPdf(const RooAbsPdf& other, const char* name=nullptr);
  ~RooAbsPdf() override;

  // Toy MC generation

  /// Keyword arguments for RooAbsPdf::generate().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct GenerateKwargs {
     RooFit::Detail::Kwarg::AllBinned AllBinned;
     RooFit::Detail::Kwarg::Asimov Asimov;
     RooFit::Detail::Kwarg::AutoBinned AutoBinned;
     RooFit::Detail::Kwarg::ExpectedData ExpectedData;
     RooFit::Detail::Kwarg::Extended Extended;
     RooFit::Detail::Kwarg::GenBinned GenBinned;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::NumEvents NumEvents;
     RooFit::Detail::Kwarg::ProtoData ProtoData;
     RooFit::Detail::Kwarg::Verbose Verbose;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&AllBinned, &Asimov, &AutoBinned, &ExpectedData, &Extended, &GenBinned, &Name, &NumEvents, &ProtoData,
                &Verbose};
        // clang-format on
     }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// See RooAbsPdf::generate(const RooArgSet&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&,const RooCmdArg&)
  /// \param[in] whatVars Set of observables to generate for each event according to this model.
  /// \param[in] nEvents How many events to generate
  /// \param arg1,arg2,arg3,arg4,arg5 Optional command arguments.
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars, Int_t nEvents, const RooCmdArg& arg1,
                       const RooCmdArg& arg2={}, const RooCmdArg& arg3={},
                       const RooCmdArg& arg4={}, const RooCmdArg& arg5={}) {
    return generate(whatVars,RooFit::NumEvents(nEvents),arg1,arg2,arg3,arg4,arg5) ;
  }
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars,
                       const RooCmdArg& arg1={},const RooCmdArg& arg2={},
                       const RooCmdArg& arg3={},const RooCmdArg& arg4={},
                       const RooCmdArg& arg5={},const RooCmdArg& arg6={}) ;
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars, GenerateKwargs const &kwargs,
                                         const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                                         const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                                         const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {})
  {
     auto a =
        RooFit::Detail::cmdArgSlots<6>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6}, "RooAbsPdf::generate()");
     return generate(whatVars, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5]);
  }
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars, Int_t nEvents, GenerateKwargs const &kwargs,
                                         const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                                         const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                                         const RooCmdArg &arg5 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<5>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5}, "RooAbsPdf::generate()");
     return generate(whatVars, nEvents, *a[0], *a[1], *a[2], *a[3], *a[4]);
  }
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars, double nEvents = 0, bool verbose=false, bool autoBinned=true,
             const char* binnedTag="", bool expectedData=false, bool extended = false) const;
  RooFit::OwningPtr<RooDataSet> generate(const RooArgSet &whatVars, const RooDataSet &prototype, Int_t nEvents= 0,
             bool verbose=false, bool randProtoOrder=false, bool resampleProto=false) const;


  class GenSpec {
  public:
    virtual ~GenSpec() ;
    GenSpec() = default;

  private:
    GenSpec(RooAbsGenContext* context, const RooArgSet& whatVars, RooDataSet* protoData, Int_t nGen, bool extended,
       bool randProto, bool resampleProto, TString dsetName, bool init=false) ;
    GenSpec(const GenSpec& other) ;

    friend class RooAbsPdf ;
    std::unique_ptr<RooAbsGenContext> _genContext;
    RooArgSet _whatVars ;
    RooDataSet* _protoData = nullptr;
    Int_t _nGen = 0;
    bool _extended = false;
    bool _randProto = false;
    bool _resampleProto = false;
    TString _dsetName ;
    bool _init = false;

    ClassDef(GenSpec,0) // Generation specification
  } ;

  /// Keyword arguments for RooAbsPdf::prepareMultiGen().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct PrepareMultiGenKwargs {
     RooFit::Detail::Kwarg::AllBinned AllBinned;
     RooFit::Detail::Kwarg::AutoBinned AutoBinned;
     RooFit::Detail::Kwarg::Extended Extended;
     RooFit::Detail::Kwarg::GenBinned GenBinned;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::NumEvents NumEvents;
     RooFit::Detail::Kwarg::ProtoData ProtoData;
     RooFit::Detail::Kwarg::Verbose Verbose;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&AllBinned, &AutoBinned, &Extended, &GenBinned, &Name, &NumEvents, &ProtoData, &Verbose};
        // clang-format on
     }
  };

  ///Prepare GenSpec configuration object for efficient generation of multiple datasets from identical specification.
  GenSpec* prepareMultiGen(const RooArgSet &whatVars,
            const RooCmdArg& arg1={},const RooCmdArg& arg2={},
            const RooCmdArg& arg3={},const RooCmdArg& arg4={},
            const RooCmdArg& arg5={},const RooCmdArg& arg6={}) ;
  GenSpec *prepareMultiGen(const RooArgSet &whatVars, PrepareMultiGenKwargs const &kwargs, const RooCmdArg &arg1 = {},
                           const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                           const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<6>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6},
                                             "RooAbsPdf::prepareMultiGen()");
     return prepareMultiGen(whatVars, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5]);
  }
  ///Generate according to GenSpec obtained from prepareMultiGen().
  RooFit::OwningPtr<RooDataSet> generate(GenSpec&) const ;

  /// Keyword arguments for RooAbsPdf::generateBinned().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct GenerateBinnedKwargs {
     RooFit::Detail::Kwarg::Asimov Asimov;
     RooFit::Detail::Kwarg::ExpectedData ExpectedData;
     RooFit::Detail::Kwarg::Extended Extended;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::NumEvents NumEvents;
     RooFit::Detail::Kwarg::Verbose Verbose;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&Asimov, &ExpectedData, &Extended, &Name, &NumEvents, &Verbose};
        // clang-format on
     }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// As RooAbsPdf::generateBinned(const RooArgSet&, const RooCmdArg&,const RooCmdArg&, const RooCmdArg&,const RooCmdArg&, const RooCmdArg&,const RooCmdArg&) const.
  /// \param[in] whatVars set
  /// \param[in] nEvents How many events to generate
  /// \param arg1,arg2,arg3,arg4,arg5 ordered arguments
  virtual RooFit::OwningPtr<RooDataHist> generateBinned(const RooArgSet &whatVars, double nEvents, const RooCmdArg& arg1,
               const RooCmdArg& arg2={}, const RooCmdArg& arg3={},
               const RooCmdArg& arg4={}, const RooCmdArg& arg5={}) const {
    return generateBinned(whatVars,RooFit::NumEvents(nEvents),arg1,arg2,arg3,arg4,arg5);
  }
  virtual RooFit::OwningPtr<RooDataHist> generateBinned(const RooArgSet &whatVars,
               const RooCmdArg& arg1={},const RooCmdArg& arg2={},
               const RooCmdArg& arg3={},const RooCmdArg& arg4={},
               const RooCmdArg& arg5={},const RooCmdArg& arg6={}) const;
  RooFit::OwningPtr<RooDataHist> generateBinned(const RooArgSet &whatVars, GenerateBinnedKwargs const &kwargs,
                                                const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                                                const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                                                const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}) const
  {
     auto a = RooFit::Detail::cmdArgSlots<6>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6},
                                             "RooAbsPdf::generateBinned()");
     return generateBinned(whatVars, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5]);
  }
  RooFit::OwningPtr<RooDataHist> generateBinned(const RooArgSet &whatVars, double nEvents,
                                                GenerateBinnedKwargs const &kwargs, const RooCmdArg &arg1 = {},
                                                const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {},
                                                const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {}) const
  {
     auto a =
        RooFit::Detail::cmdArgSlots<5>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5}, "RooAbsPdf::generateBinned()");
     return generateBinned(whatVars, nEvents, *a[0], *a[1], *a[2], *a[3], *a[4]);
  }
  virtual RooFit::OwningPtr<RooDataHist> generateBinned(const RooArgSet &whatVars, double nEvents, bool expectedData=false, bool extended=false) const;

  virtual RooFit::OwningPtr<RooDataSet> generateSimGlobal(const RooArgSet& whatVars, Int_t nEvents) ;

  /// Keyword arguments for RooAbsPdf::plotOn().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct PlotOnKwargs {
     RooFit::Detail::Kwarg::AddTo AddTo;
     RooFit::Detail::Kwarg::Asymmetry Asymmetry;
     RooFit::Detail::Kwarg::Components Components;
     RooFit::Detail::Kwarg::DrawOption DrawOption;
     RooFit::Detail::Kwarg::EvalErrorValue EvalErrorValue;
     RooFit::Detail::Kwarg::FillColor FillColor;
     RooFit::Detail::Kwarg::FillStyle FillStyle;
     RooFit::Detail::Kwarg::Invisible Invisible;
     RooFit::Detail::Kwarg::LineColor LineColor;
     RooFit::Detail::Kwarg::LineStyle LineStyle;
     RooFit::Detail::Kwarg::LineWidth LineWidth;
     RooFit::Detail::Kwarg::MarkerColor MarkerColor;
     RooFit::Detail::Kwarg::MarkerSize MarkerSize;
     RooFit::Detail::Kwarg::MarkerStyle MarkerStyle;
     RooFit::Detail::Kwarg::MoveToBack MoveToBack;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::NormRange NormRange;
     RooFit::Detail::Kwarg::Normalization Normalization;
     RooFit::Detail::Kwarg::NumCPU NumCPU;
     RooFit::Detail::Kwarg::Precision Precision;
     RooFit::Detail::Kwarg::PrintEvalErrors PrintEvalErrors;
     RooFit::Detail::Kwarg::ProjWData ProjWData;
     RooFit::Detail::Kwarg::Project Project;
     RooFit::Detail::Kwarg::ProjectionRange ProjectionRange;
     RooFit::Detail::Kwarg::Range Range;
     RooFit::Detail::Kwarg::ShiftToZero ShiftToZero;
     RooFit::Detail::Kwarg::ShowProgress ShowProgress;
     RooFit::Detail::Kwarg::Slice Slice;
     RooFit::Detail::Kwarg::VLines VLines;
     RooFit::Detail::Kwarg::VisualizeError VisualizeError;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&AddTo, &Asymmetry, &Components, &DrawOption, &EvalErrorValue, &FillColor, &FillStyle, &Invisible,
                &LineColor, &LineStyle, &LineWidth, &MarkerColor, &MarkerSize, &MarkerStyle, &MoveToBack, &Name,
                &NormRange, &Normalization, &NumCPU, &Precision, &PrintEvalErrors, &ProjWData, &Project,
                &ProjectionRange, &Range, &ShiftToZero, &ShowProgress, &Slice, &VLines, &VisualizeError};
        // clang-format on
     }
  };

  ///Helper calling plotOn(RooPlot*, RooLinkedList&) const
  RooPlot* plotOn(RooPlot* frame,
           const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
           const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
           const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
           const RooCmdArg& arg7={}, const RooCmdArg& arg8={},
           const RooCmdArg& arg9={}, const RooCmdArg& arg10={}
              ) const override {
    return RooAbsReal::plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10) ;
  }
  RooPlot* plotOn(RooPlot* frame, RooLinkedList& cmdList) const override ;
  RooPlot *plotOn(RooPlot *frame, PlotOnKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                  const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                  const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {},
                  const RooCmdArg &arg9 = {}, const RooCmdArg &arg10 = {}) const
  {
     return plotOn(frame, *RooFit::Detail::createCmdList(
                             kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8, &arg9, &arg10}));
  }

  /// Keyword arguments for RooAbsPdf::paramOn().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct ParamOnKwargs {
     RooFit::Detail::Kwarg::Format Format;
     RooFit::Detail::Kwarg::Label Label;
     RooFit::Detail::Kwarg::Layout Layout;
     RooFit::Detail::Kwarg::Parameters Parameters;
     RooFit::Detail::Kwarg::ShowConstants ShowConstants;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&Format, &Label, &Layout, &Parameters, &ShowConstants};
        // clang-format on
     }
  };

  /// Add a box with parameter values (and errors) to the specified frame
  virtual RooPlot* paramOn(RooPlot* frame,
                           const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
                           const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
                           const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
                           const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooPlot *paramOn(RooPlot *frame, ParamOnKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                   const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                   const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooAbsPdf::paramOn()");
     return paramOn(frame, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }

  // Built-in generator support
  virtual Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, bool staticInitOK=true) const;
  virtual void initGenerator(Int_t code) ;
  virtual void generateEvent(Int_t code);
  virtual bool isDirectGenSafe(const RooAbsArg& arg) const ;

  // Configuration of MC generators used for this pdf
  const RooNumGenConfig* getGeneratorConfig() const ;
  static RooNumGenConfig* defaultGeneratorConfig()  ;
  RooNumGenConfig* specialGeneratorConfig() const ;
  RooNumGenConfig* specialGeneratorConfig(bool createOnTheFly) ;
  void setGeneratorConfig() ;
  void setGeneratorConfig(const RooNumGenConfig& config) ;

  template <typename... CmdArgs_t>
  RooFit::OwningPtr<RooFitResult> fitTo(RooAbsData& data, CmdArgs_t const&... cmdArgs)
  {
    return RooFit::makeOwningPtr(fitToImpl(data, *RooFit::Detail::createCmdList(&cmdArgs...)));
  }

  RooFit::OwningPtr<RooFitResult>
  fitTo(RooAbsData &data, RooFit::FitToKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
        const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {},
        const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {});

  /// Keyword arguments for RooAbsPdf::createNLL().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct CreateNLLKwargs {
     RooFit::Detail::Kwarg::CloneData CloneData;
     RooFit::Detail::Kwarg::ConditionalObservables ConditionalObservables;
     RooFit::Detail::Kwarg::Constrain Constrain;
     RooFit::Detail::Kwarg::EvalBackend EvalBackend;
     RooFit::Detail::Kwarg::Extended Extended;
     RooFit::Detail::Kwarg::ExternalConstraints ExternalConstraints;
     RooFit::Detail::Kwarg::GlobalObservables GlobalObservables;
     RooFit::Detail::Kwarg::GlobalObservablesSource GlobalObservablesSource;
     RooFit::Detail::Kwarg::GlobalObservablesTag GlobalObservablesTag;
     RooFit::Detail::Kwarg::IntegrateBins IntegrateBins;
     RooFit::Detail::Kwarg::ModularL ModularL;
     RooFit::Detail::Kwarg::NumCPU NumCPU;
     RooFit::Detail::Kwarg::Offset Offset;
     RooFit::Detail::Kwarg::ProjectedObservables ProjectedObservables;
     RooFit::Detail::Kwarg::Range Range;
     RooFit::Detail::Kwarg::SplitRange SplitRange;
     RooFit::Detail::Kwarg::SumCoefRange SumCoefRange;
     RooFit::Detail::Kwarg::Verbose Verbose;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&CloneData, &ConditionalObservables, &Constrain, &EvalBackend, &Extended, &ExternalConstraints,
                &GlobalObservables, &GlobalObservablesSource, &GlobalObservablesTag, &IntegrateBins, &ModularL,
                &NumCPU, &Offset, &ProjectedObservables, &Range, &SplitRange, &SumCoefRange, &Verbose};
        // clang-format on
     }
  };

  template <typename... CmdArgs_t>
  RooFit::OwningPtr<RooAbsReal> createNLL(RooAbsData& data, CmdArgs_t const&... cmdArgs)
  {
    return RooFit::makeOwningPtr(createNLLImpl(data, *RooFit::Detail::createCmdList(&cmdArgs...)));
  }
  RooFit::OwningPtr<RooAbsReal>
  createNLL(RooAbsData &data, CreateNLLKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
            const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
            const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
  {
     return RooFit::makeOwningPtr(createNLLImpl(
        data, *RooFit::Detail::createCmdList(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8})));
  }

  // Constraint management. Interface to retrieve constraint terms on this pdf. Default implementation returns null.
  virtual std::unique_ptr<RooArgSet> getConstraints(const RooArgSet & /*observables*/,
                                                    RooArgSet const & /*constrainedParams*/,
                                                    RooArgSet & /*pdfParams*/) const
  {
     return nullptr;
  }

  std::unique_ptr<RooArgSet>
  getAllConstraints(const RooArgSet &observables, RooArgSet &constrainedParams, bool stripDisconnected = true) const;

  // Project p.d.f into lower dimensional p.d.f
  virtual RooAbsPdf* createProjection(const RooArgSet& iset) ;

  virtual double getCorrection() const;

  // Create cumulative density function from p.d.f
  /// Keyword arguments for RooAbsPdf::createCdf().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct CreateCdfKwargs {
     RooFit::Detail::Kwarg::ScanAllCdf ScanAllCdf;
     RooFit::Detail::Kwarg::ScanNoCdf ScanNoCdf;
     RooFit::Detail::Kwarg::ScanNumCdf ScanNumCdf;
     RooFit::Detail::Kwarg::ScanParameters ScanParameters;
     RooFit::Detail::Kwarg::SupNormSet SupNormSet;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&ScanAllCdf, &ScanNoCdf, &ScanNumCdf, &ScanParameters, &SupNormSet};
        // clang-format on
     }
  };

  RooFit::OwningPtr<RooAbsReal> createCdf(const RooArgSet& iset, const RooArgSet& nset=RooArgSet()) ;
  RooFit::OwningPtr<RooAbsReal> createCdf(const RooArgSet& iset, const RooCmdArg& arg1, const RooCmdArg& arg2={},
         const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
         const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
         const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) ;
  RooFit::OwningPtr<RooAbsReal> createCdf(const RooArgSet &iset, CreateCdfKwargs const &kwargs,
                                          const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                                          const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                                          const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {},
                                          const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooAbsPdf::createCdf()");
     return createCdf(iset, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }
  RooFit::OwningPtr<RooAbsReal> createScanCdf(const RooArgSet& iset, const RooArgSet& nset, Int_t numScanBins, Int_t intOrder) ;

  // Function evaluation support
  double getValV(const RooArgSet* set=nullptr) const override ;
  virtual double getLogVal(const RooArgSet* set=nullptr) const ;

  void getLogProbabilities(std::span<const double> pdfValues, double * output) const;

  /// \copydoc getNorm(const RooArgSet*) const
  double getNorm(const RooArgSet& nset) const {
    return getNorm(&nset) ;
  }
  virtual double getNorm(const RooArgSet* set=nullptr) const ;

  virtual void resetErrorCounters(Int_t resetValue=10) ;
  void setTraceCounter(Int_t value, bool allNodes=false) ;

  double analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=nullptr) const override ;

  /// Shows if a PDF is self-normalized, which means that no attempt is made to add a normalization term.
  /// Always returns false, unless a PDF overrides this function.
  virtual bool selfNormalized() const {
    return false ;
  }

  // Support for extended maximum likelihood, switched off by default
  enum ExtendMode { CanNotBeExtended, CanBeExtended, MustBeExtended } ;
  /// Returns ability of PDF to provide extended likelihood terms. Possible
  /// answers are in the enumerator RooAbsPdf::ExtendMode.
  /// This default implementation always returns CanNotBeExtended.
  virtual ExtendMode extendMode() const { return CanNotBeExtended; }
  /// If true, PDF can provide extended likelihood term.
  inline bool canBeExtended() const {
    return (extendMode() != CanNotBeExtended) ;
  }
  /// If true PDF must provide extended likelihood term.
  inline bool mustBeExtended() const {
    return (extendMode() == MustBeExtended) ;
  }
  /// Return expected number of events to be used in calculation of extended
  /// likelihood.
  virtual double expectedEvents(const RooArgSet* nset) const ;
  /// Return expected number of events to be used in calculation of extended
  /// likelihood. This function should not be overridden, as it just redirects
  /// to the actual virtual function but takes a RooArgSet reference instead of
  /// pointer. \see expectedEvents(const RooArgSet*) const
  double expectedEvents(const RooArgSet& nset) const {
    return expectedEvents(&nset) ;
  }

  virtual std::unique_ptr<RooAbsReal> createExpectedEventsFunc(const RooArgSet* nset) const;

  // Printing interface (human readable)
  void printValue(std::ostream& os) const override ;
  void printMultiline(std::ostream& os, Int_t contents, bool verbose=false, TString indent="") const override ;

  static void verboseEval(Int_t stat) ;
  static int verboseEval() ;

  double extendedTerm(double sumEntries, double expected, double sumEntriesW2=0.0, bool doOffset=false) const;
  double extendedTerm(double sumEntries, RooArgSet const* nset, double sumEntriesW2=0.0, bool doOffset=false) const;
  double extendedTerm(RooAbsData const& data, bool weightSquared, bool doOffset=false) const;

  void setNormRange(const char* rangeName) ;
  const char* normRange() const {
    return _normRange.Length()>0 ? _normRange.Data() : nullptr ;
  }
  void setNormRangeOverride(const char* rangeName) ;

  const RooAbsReal* getNormIntegral(const RooArgSet& nset) const { return getNormObj(nullptr,&nset,nullptr) ; }

  virtual const RooAbsReal* getNormObj(const RooArgSet* set, const RooArgSet* iset, const TNamed* rangeName=nullptr) const ;

  virtual RooAbsGenContext* binnedGenContext(const RooArgSet &vars, bool verbose= false) const ;

  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=nullptr,
                                  const RooArgSet* auxProto=nullptr, bool verbose= false) const ;

  virtual RooAbsGenContext* autoGenContext(const RooArgSet &vars, const RooDataSet* prototype=nullptr, const RooArgSet* auxProto=nullptr,
                  bool verbose=false, bool autoBinned=true, const char* binnedTag="") const ;

  std::unique_ptr<RooAbsArg> compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext & ctx) const override;

private:

  std::unique_ptr<RooDataSet> generate(RooAbsGenContext& context, const RooArgSet& whatVars, const RooDataSet* prototype,
             double nEvents, bool verbose, bool randProtoOrder, bool resampleProto, bool skipInit=false,
             bool extended=false) const ;

  // Implementation version
  virtual RooPlot* paramOn(RooPlot* frame, const RooArgSet& params, bool showConstants=false,
                           const char *label= "", double xmin=0.65,
                           double xmax= 0.99,double ymax=0.95, const RooCmdArg* formatCmd=nullptr) ;

  void logBatchComputationErrors(std::span<const double>& outputs, std::size_t begin) const;
  bool traceEvalPdf(double value) const;

  /// Setter for the _normSet member, which should never be set directly.
  inline void setActiveNormSet(RooArgSet const* normSet) const {
    _normSet = normSet;
    // Also store the unique ID of the _normSet. This makes it possible to
    // detect if the pointer was invalidated.
    _normSetId = RooFit::getUniqueId(normSet);
  }

protected:

  virtual std::unique_ptr<RooAbsReal> createNLLImpl(RooAbsData& data, const RooLinkedList& cmdList);
  virtual std::unique_ptr<RooFitResult> fitToImpl(RooAbsData& data, const RooLinkedList& cmdList);

  /// Checks if `normSet` is the currently active normalization set of this
  /// PDF, meaning is exactly the same object as the one the `_normSet` member
  /// points to (or both are `nullptr`).
  inline bool isActiveNormSet(RooArgSet const* normSet) const {
    return RooFit::getUniqueId(normSet).value() == _normSetId;
  }

  RooPlot *plotOn(RooPlot *frame, PlotOpt o) const override;

  friend class RooMCStudy ;

  Int_t* randomizeProtoOrder(Int_t nProto,Int_t nGen,bool resample=false) const ;

  // This also forces the definition of a copy ctor in derived classes
  RooAbsPdf(const RooAbsPdf& other, const char* name = nullptr);

  static Int_t _verboseEval ;

  virtual bool syncNormalization(const RooArgSet* dset, bool adjustProxies=true) const ;

  double _rawValue = 0;
  mutable RooAbsReal* _norm = nullptr; ///<! Normalization integral (owned by _normMgr)
  mutable RooArgSet const* _normSet = nullptr; ///<! Normalization set with for above integral

  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem(RooAbsReal& norm) : _norm(&norm) {} ;
    ~CacheElem() override ;
    RooArgList containedArgs(Action) override { return RooArgList(*_norm) ; }
    std::unique_ptr<RooAbsReal> _norm;
  } ;
  mutable RooObjCacheManager _normMgr ; ///<! The cache manager

  bool redirectServersHook(const RooAbsCollection & newServerList, bool mustReplaceAll,
                                   bool nameChange, bool isRecursiveStep) override;

  mutable Int_t _errorCount = 0; ///< Number of errors remaining to print
  Int_t _traceCount = 0;         ///< Number of traces remaining to print
  Int_t _negCount = 0;           ///< Number of negative probabilities remaining to print

  bool _selectComp = false; ///< Component selection flag for RooAbsPdf::plotCompOn

  std::unique_ptr<RooNumGenConfig> _specGeneratorConfig ; ///<! MC generator configuration specific for this object

  TString _normRange ; ///< Normalization range
  static TString _normRangeOverride ;

private:
  mutable RooFit::UniqueId<RooArgSet>::Value_t _normSetId = RooFit::UniqueId<RooArgSet>::nullval; ///<! Unique ID of the currently-active normalization set

  friend class RooAbsReal;

  ClassDefOverride(RooAbsPdf,5) // Abstract PDF with normalization support
};




#endif
