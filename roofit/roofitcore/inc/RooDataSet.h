/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDataSet.h,v 1.59 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_DATA_SET
#define ROO_DATA_SET

class RooAbsRealLValue;
class RooCategory;
class RooDataHist;
class RooRealVar;

class TDirectory;

#include "RooAbsData.h"
#include "RooDirItem.h"
#include "RooFit/Kwargs.h"

#include <list>
#include <string_view>

class RooDataSet : public RooAbsData, public RooDirItem {
public:

  // Constructors, factory methods etc.
  RooDataSet() ;

  // Universal constructor
  /// Keyword arguments for the RooDataSet constructor.
  /// Each member corresponds to the RooFit command argument of the same name.
  struct Kwargs {
     RooFit::Detail::Kwarg::Cut Cut;
     RooFit::Detail::Kwarg::CutRange CutRange;
     RooFit::Detail::Kwarg::GlobalObservables GlobalObservables;
     RooFit::Detail::Kwarg::Import Import;
     RooFit::Detail::Kwarg::ImportFromFile ImportFromFile;
     RooFit::Detail::Kwarg::Index Index;
     RooFit::Detail::Kwarg::Link Link;
     RooFit::Detail::Kwarg::OwnLinked OwnLinked;
     RooFit::Detail::Kwarg::StoreAsymError StoreAsymError;
     RooFit::Detail::Kwarg::StoreError StoreError;
     RooFit::Detail::Kwarg::WeightVar WeightVar;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&Cut, &CutRange, &GlobalObservables, &Import, &ImportFromFile, &Index, &Link, &OwnLinked,
                &StoreAsymError, &StoreError, &WeightVar};
        // clang-format on
     }
  };

  RooDataSet(RooStringView name, RooStringView title, const RooArgSet& vars, const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
             const RooCmdArg& arg3={}, const RooCmdArg& arg4={},const RooCmdArg& arg5={},
             const RooCmdArg& arg6={},const RooCmdArg& arg7={},const RooCmdArg& arg8={}) ;
  RooDataSet(RooStringView name, RooStringView title, const RooArgSet &vars, Kwargs const &kwargs,
             const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {},
             const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {},
             const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
     : RooDataSet(RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                                 "RooDataSet::RooDataSet()"),
                  name, title, vars)
  {
  }

  RooDataSet(RooDataSet const & other, const char* newname=nullptr) ;
  TObject* Clone(const char* newname = "") const override {
    return new RooDataSet(*this, newname && newname[0] != '\0' ? newname : GetName());
  }
  ~RooDataSet() override ;

  RooFit::OwningPtr<RooAbsData> emptyClone(const char* newName=nullptr, const char* newTitle=nullptr, const RooArgSet* vars=nullptr, const char* wgtVarName=nullptr) const override;

  RooFit::OwningPtr<RooDataHist> binnedClone(const char* newName=nullptr, const char* newTitle=nullptr) const ;

  double sumEntries() const override;
  double sumEntries(const char* cutSpec, const char* cutRange=nullptr) const override;

  /// Keyword arguments for RooDataSet::plotOnXY().
  /// Each member corresponds to the RooFit command argument of the same name.
  struct PlotOnXYKwargs {
     RooFit::Detail::Kwarg::DrawOption DrawOption;
     RooFit::Detail::Kwarg::FillColor FillColor;
     RooFit::Detail::Kwarg::FillStyle FillStyle;
     RooFit::Detail::Kwarg::Invisible Invisible;
     RooFit::Detail::Kwarg::LineColor LineColor;
     RooFit::Detail::Kwarg::LineStyle LineStyle;
     RooFit::Detail::Kwarg::LineWidth LineWidth;
     RooFit::Detail::Kwarg::MarkerColor MarkerColor;
     RooFit::Detail::Kwarg::MarkerSize MarkerSize;
     RooFit::Detail::Kwarg::MarkerStyle MarkerStyle;
     RooFit::Detail::Kwarg::Name Name;
     RooFit::Detail::Kwarg::Rescale Rescale;
     RooFit::Detail::Kwarg::YVar YVar;

     std::vector<RooCmdArg const *> cmdArgs() const
     {
        // clang-format off
        return {&DrawOption, &FillColor, &FillStyle, &Invisible, &LineColor, &LineStyle, &LineWidth, &MarkerColor,
                &MarkerSize, &MarkerStyle, &Name, &Rescale, &YVar};
        // clang-format on
     }
  };

  virtual RooPlot* plotOnXY(RooPlot* frame,
             const RooCmdArg& arg1={}, const RooCmdArg& arg2={},
             const RooCmdArg& arg3={}, const RooCmdArg& arg4={},
             const RooCmdArg& arg5={}, const RooCmdArg& arg6={},
             const RooCmdArg& arg7={}, const RooCmdArg& arg8={}) const ;
  RooPlot *plotOnXY(RooPlot *frame, PlotOnXYKwargs const &kwargs, const RooCmdArg &arg1 = {},
                    const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
                    const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {},
                    const RooCmdArg &arg8 = {}) const
  {
     auto a = RooFit::Detail::cmdArgSlots<8>(kwargs, {&arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8},
                                             "RooDataSet::plotOnXY()");
     return plotOnXY(frame, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7]);
  }

  /// Read data from a text file and create a dataset from it.
  /// The possible options are: (D)ebug, (Q)uiet.
  static RooDataSet *read(const char *filename, const RooArgList &variables,
           const char *opts= "", const char* commonPath="",
           const char *indexCatName=nullptr) ;
  bool write(const char* filename) const;
  bool write(std::ostream & ofs) const;


  bool isWeighted() const override;
  bool isNonPoissonWeighted() const override;

  double weight() const override;
  /// Returns a pointer to the weight variable (if set).
  RooRealVar* weightVar() const { return _wgtVar; }
  double weightSquared() const override;
  void weightError(double& lo, double& hi,ErrorType etype=SumW2) const override;
  double weightError(ErrorType etype=SumW2) const override;

  const RooArgSet* get(Int_t index) const override;
  const RooArgSet* get() const override;

  std::span<const double> getWeightBatch(std::size_t first, std::size_t len, bool sumW2) const override;

  /// Add one ore more rows of data
  void add(const RooArgSet& row, double weight, double weightError);
  void add(const RooArgSet& row, double weight=1.0) override { add(row, weight, 0.0); }
  virtual void add(const RooArgSet& row, double weight, double weightErrorLo, double weightErrorHi);

  virtual void addFast(const RooArgSet& row, double weight=1.0, double weightError=0.0);

  void append(RooDataSet& data) ;
  bool merge(RooDataSet* data1, RooDataSet* data2=nullptr, RooDataSet* data3=nullptr,
           RooDataSet* data4=nullptr, RooDataSet* data5=nullptr, RooDataSet* data6=nullptr) ;
  bool merge(std::list<RooDataSet*> dsetList) ;

  virtual RooAbsArg* addColumn(RooAbsArg& var, bool adjustRange=true) ;

  void printMultiline(std::ostream& os, Int_t contents, bool verbose=false, TString indent="") const override;
  void printArgs(std::ostream& os) const override;
  void printValue(std::ostream& os) const override;

  /// Print the contents of the dataset to the specified output stream.
  void printContents(std::ostream& os = std::cout) const override;

  void SetName(const char *name) override;
  void SetNameTitle(const char *name, const char* title) override;

  static void cleanup();

  void convertToTreeStore() override;

protected:

  friend class RooProdGenContext ;

  void initialize(const char* wgtVarName) ;

  // Cache copy feature is not publicly accessible
  std::unique_ptr<RooAbsData> reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange=nullptr,
                        std::size_t nStart=0, std::size_t nStop = std::numeric_limits<std::size_t>::max()) const override;

  RooArgSet _varsNoWgt;          ///< Vars without weight variable
  RooRealVar *_wgtVar = nullptr; ///< Pointer to weight variable (if set)

private:
   // Delegated to by the constructor that takes kwargs.
   RooDataSet(std::array<RooCmdArg const *, 8> const &a, RooStringView name, RooStringView title, const RooArgSet &vars)
      : RooDataSet(name, title, vars, *a[0], *a[1], *a[2], *a[3], *a[4], *a[5], *a[6], *a[7])
   {
   }

  void loadValuesFromSlices(RooCategory &indexCat, std::map<std::string, RooAbsData *> const &slices,
                            const char *rangeName, RooFormulaVar const *cutVar, const char *cutSpec);

  unsigned short _errorMsgCount{0}; ///<! Counter to silence error messages when filling dataset.
  bool _doWeightErrorCheck{true};   ///<! When adding events with weights, check that weights can actually be stored.

  mutable std::unique_ptr<std::vector<double>> _sumW2Buffer; ///<! Buffer for sumW2 in case a batch of values is requested.

  ClassDefOverride(RooDataSet,2) // Unbinned data set
};

#endif
