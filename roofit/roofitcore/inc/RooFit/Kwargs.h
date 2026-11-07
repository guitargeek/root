/*
 * Project: RooFit
 * Authors:
 *   Jonas Rembser, CERN 2026
 *
 * Copyright (c) 2026, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_Kwargs_h
#define RooFit_Kwargs_h

#include <RooCmdArg.h>
#include <RooGlobalFunc.h>
#include <RooLinkedList.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
\file Kwargs.h
\brief Keyword arguments structs, the alternative to RooFit command arguments.

Many RooFit functions are configured with RooFit command arguments, like
RooAbsPdf::fitTo(). Most of them also have an overload that takes a struct with
a data member for each command argument they accept, named `kwargs` because it
corresponds to Python keyword arguments:

~~~{.py}
pdf.fitTo(data, Range="sideband", Save=True, PrintLevel=-1)
~~~

In C++20, designated initializers give the same with the members in the order
of declaration, which is alphabetical:

~~~{.cpp}
pdf.fitTo(data, {.PrintLevel = -1, .Range = "sideband", .Save = true});
~~~

The members can be initialized from the same arguments that the RooFit command
argument function of the same name takes, e.g. `.Range = {0., 10.}` is
equivalent to `RooFit::Range(0., 10.)`. Members that are not initialized are
ignored. Command arguments can be passed in addition, after the struct.
*/

namespace RooFit {
namespace Detail {

/// Collect the non-empty command arguments, first the given ones and then the
/// members of a kwargs struct.
template <class Kwargs_t>
std::vector<RooCmdArg const *> collectCmdArgs(Kwargs_t const &kwargs, std::initializer_list<RooCmdArg const *> args)
{
   std::vector<RooCmdArg const *> out;
   for (RooCmdArg const *arg : args) {
      if (arg->opcode())
         out.push_back(arg);
   }
   for (RooCmdArg const *arg : kwargs.cmdArgs()) {
      if (arg->opcode())
         out.push_back(arg);
   }
   return out;
}

/// Put the command arguments from collectCmdArgs() in a RooLinkedList.
template <class Kwargs_t>
std::unique_ptr<RooLinkedList> createCmdList(Kwargs_t const &kwargs, std::initializer_list<RooCmdArg const *> args)
{
   auto cmdList = std::make_unique<RooLinkedList>();
   for (RooCmdArg const *arg : collectCmdArgs(kwargs, args)) {
      cmdList->Add(const_cast<RooCmdArg *>(arg));
   }
   return cmdList;
}

/// Distribute the command arguments from collectCmdArgs() over the N RooCmdArg
/// parameters of a function that doesn't take a RooLinkedList.
template <std::size_t N, class Kwargs_t>
std::array<RooCmdArg const *, N>
cmdArgSlots(Kwargs_t const &kwargs, std::initializer_list<RooCmdArg const *> args, const char *funcName)
{
   std::vector<RooCmdArg const *> cmdArgs = collectCmdArgs(kwargs, args);
   if (cmdArgs.size() > N) {
      throw std::invalid_argument(std::string{funcName} + " takes at most " + std::to_string(N) +
                                  " command arguments, but got " + std::to_string(cmdArgs.size()));
   }
   std::array<RooCmdArg const *, N> out;
   out.fill(&RooCmdArg::none());
   std::copy(cmdArgs.begin(), cmdArgs.end(), out.begin());
   return out;
}

/// Types of the data members of kwargs structs, one for each RooFit command
/// argument. They are command arguments themselves, with constructors that
/// mirror the overloads of the function with the same name. A default
/// constructed instance is an empty command argument, which is ignored.
namespace Kwarg {

/// \see RooFit::AddTo()
struct AddTo : RooCmdArg {
   AddTo() = default;
   AddTo(const char *name, double wgtSel = 1.0, double wgtOther = 1.0)
      : RooCmdArg{RooFit::AddTo(name, wgtSel, wgtOther)}
   {
   }
};

/// \see RooFit::AllBinned()
struct AllBinned : RooCmdArg {
   AllBinned() = default;
   AllBinned(bool flag) : RooCmdArg{flag ? RooFit::AllBinned() : RooCmdArg{}} {}
};

/// \see RooFit::Asimov()
struct Asimov : RooCmdArg {
   Asimov() = default;
   Asimov(bool flag) : RooCmdArg{RooFit::Asimov(flag)} {}
};

/// \see RooFit::Asymmetry()
struct Asymmetry : RooCmdArg {
   Asymmetry() = default;
   Asymmetry(const RooCategory &cat) : RooCmdArg{RooFit::Asymmetry(cat)} {}
};

/// \see RooFit::AsymptoticError()
struct AsymptoticError : RooCmdArg {
   AsymptoticError() = default;
   AsymptoticError(bool flag) : RooCmdArg{RooFit::AsymptoticError(flag)} {}
};

/// \see RooFit::AutoBinned()
struct AutoBinned : RooCmdArg {
   AutoBinned() = default;
   AutoBinned(bool flag) : RooCmdArg{RooFit::AutoBinned(flag)} {}
};

/// \see RooFit::AutoBinning()
struct AutoBinning : RooCmdArg {
   AutoBinning() = default;
   AutoBinning(Int_t nbins, double marginFactor = 0.1) : RooCmdArg{RooFit::AutoBinning(nbins, marginFactor)} {}
};

/// \see RooFit::AutoPrecision()
struct AutoPrecision : RooCmdArg {
   AutoPrecision() = default;
   AutoPrecision(Int_t ndigit) : RooCmdArg{RooFit::AutoPrecision(ndigit)} {}
};

/// \see RooFit::AutoRange()
struct AutoRange : RooCmdArg {
   AutoRange() = default;
   AutoRange(const RooAbsData &data, double marginFactor = 0.1) : RooCmdArg{RooFit::AutoRange(data, marginFactor)} {}
};

/// \see RooFit::AutoSymBinning()
struct AutoSymBinning : RooCmdArg {
   AutoSymBinning() = default;
   AutoSymBinning(Int_t nbins, double marginFactor = 0.1) : RooCmdArg{RooFit::AutoSymBinning(nbins, marginFactor)} {}
};

/// \see RooFit::AutoSymRange()
struct AutoSymRange : RooCmdArg {
   AutoSymRange() = default;
   AutoSymRange(const RooAbsData &data, double marginFactor = 0.1) : RooCmdArg{RooFit::AutoSymRange(data, marginFactor)}
   {
   }
};

/// \see RooFit::AxisLabel()
struct AxisLabel : RooCmdArg {
   AxisLabel() = default;
   AxisLabel(const char *name) : RooCmdArg{RooFit::AxisLabel(name)} {}
};

/// \see RooFit::BaseClassName()
struct BaseClassName : RooCmdArg {
   BaseClassName() = default;
   BaseClassName(const char *name) : RooCmdArg{RooFit::BaseClassName(name)} {}
};

/// \see RooFit::Binned()
struct Binned : RooCmdArg {
   Binned() = default;
   Binned(bool flag) : RooCmdArg{RooFit::Binned(flag)} {}
};

/// \see RooFit::Binning()
struct Binning : RooCmdArg {
   Binning() = default;
   Binning(const RooAbsBinning &binning) : RooCmdArg{RooFit::Binning(binning)} {}
   Binning(const char *binningName) : RooCmdArg{RooFit::Binning(binningName)} {}
   Binning(int nBins, double xlo = 0.0, double xhi = 0.0) : RooCmdArg{RooFit::Binning(nBins, xlo, xhi)} {}
};

/// \see RooFit::Bins()
struct Bins : RooCmdArg {
   Bins() = default;
   Bins(Int_t nbin) : RooCmdArg{RooFit::Bins(nbin)} {}
};

/// \see RooFit::ClassName()
struct ClassName : RooCmdArg {
   ClassName() = default;
   ClassName(const char *name) : RooCmdArg{RooFit::ClassName(name)} {}
};

/// \see RooFit::CloneData()
struct CloneData : RooCmdArg {
   CloneData() = default;
   CloneData(bool flag) : RooCmdArg{RooFit::CloneData(flag)} {}
};

/// \see RooFit::Color()
struct Color : RooCmdArg {
   Color() = default;
   Color(TColorNumber color) : RooCmdArg{RooFit::Color(color)} {}
};

/// \see RooFit::Columns()
struct Columns : RooCmdArg {
   Columns() = default;
   Columns(Int_t ncol) : RooCmdArg{RooFit::Columns(ncol)} {}
};

/// \see RooFit::Components()
struct Components : RooCmdArg {
   Components() = default;
   Components(const char *compSpec) : RooCmdArg{RooFit::Components(compSpec)} {}
   Components(const RooArgSet &args) : RooCmdArg{RooFit::Components(args)} {}
};

/// \see RooFit::Conditional()
struct Conditional : RooCmdArg {
   Conditional() = default;
   Conditional(const RooArgSet &pdfSet, const RooArgSet &depSet, bool depsAreCond = false)
      : RooCmdArg{RooFit::Conditional(pdfSet, depSet, depsAreCond)}
   {
   }
};

/// \see RooFit::ConditionalObservables()
struct ConditionalObservables : RooCmdArg {
   ConditionalObservables() = default;
   ConditionalObservables(const RooArgSet &args) : RooCmdArg{RooFit::ConditionalObservables(args)} {}
};

/// \see RooFit::Constrain()
struct Constrain : RooCmdArg {
   Constrain() = default;
   Constrain(const RooArgSet &params) : RooCmdArg{RooFit::Constrain(params)} {}
};

/// \see RooFit::Cut()
struct Cut : RooCmdArg {
   Cut() = default;
   Cut(const char *cutSpec) : RooCmdArg{RooFit::Cut(cutSpec)} {}
   Cut(const RooFormulaVar &cutVar) : RooCmdArg{RooFit::Cut(cutVar)} {}
};

/// \see RooFit::CutRange()
struct CutRange : RooCmdArg {
   CutRange() = default;
   CutRange(const char *rangeName) : RooCmdArg{RooFit::CutRange(rangeName)} {}
};

/// \see RooFit::DataError()
struct DataError : RooCmdArg {
   DataError() = default;
   DataError(Int_t etype) : RooCmdArg{RooFit::DataError(etype)} {}
   DataError(std::string const &etype) : RooCmdArg{RooFit::DataError(etype)} {}
   DataError(std::nullptr_t) : RooCmdArg{RooFit::DataError("None")} {}
};

/// \see RooFit::DrawOption()
struct DrawOption : RooCmdArg {
   DrawOption() = default;
   DrawOption(const char *opt) : RooCmdArg{RooFit::DrawOption(opt)} {}
};

/// \see RooFit::Efficiency()
struct Efficiency : RooCmdArg {
   Efficiency() = default;
   Efficiency(const RooCategory &cat) : RooCmdArg{RooFit::Efficiency(cat)} {}
};

/// \see RooFit::Embedded()
struct Embedded : RooCmdArg {
   Embedded() = default;
   Embedded(bool flag) : RooCmdArg{RooFit::Embedded(flag)} {}
};

/// \see RooFit::EvalBackend()
struct EvalBackend : RooCmdArg {
   EvalBackend() = default;
   EvalBackend(RooFit::EvalBackend const &backend) : RooCmdArg{backend} {}
   EvalBackend(RooFit::EvalBackend::Value value) : RooCmdArg{RooFit::EvalBackend(value)} {}
   EvalBackend(std::string const &name) : RooCmdArg{RooFit::EvalBackend(name)} {}
};

/// \see RooFit::EvalErrorValue()
struct EvalErrorValue : RooCmdArg {
   EvalErrorValue() = default;
   EvalErrorValue(double value) : RooCmdArg{RooFit::EvalErrorValue(value)} {}
};

/// \see RooFit::EvalErrorWall()
struct EvalErrorWall : RooCmdArg {
   EvalErrorWall() = default;
   EvalErrorWall(bool flag) : RooCmdArg{RooFit::EvalErrorWall(flag)} {}
};

/// \see RooFit::EventRange()
struct EventRange : RooCmdArg {
   EventRange() = default;
   EventRange(Int_t nStart, Int_t nStop) : RooCmdArg{RooFit::EventRange(nStart, nStop)} {}
};

/// \see RooFit::ExpectedData()
struct ExpectedData : RooCmdArg {
   ExpectedData() = default;
   ExpectedData(bool flag) : RooCmdArg{RooFit::ExpectedData(flag)} {}
};

/// \see RooFit::Extended()
struct Extended : RooCmdArg {
   Extended() = default;
   Extended(bool flag) : RooCmdArg{RooFit::Extended(flag)} {}
};

/// \see RooFit::ExternalConstraints()
struct ExternalConstraints : RooCmdArg {
   ExternalConstraints() = default;
   ExternalConstraints(const RooArgSet &constraintPdfs) : RooCmdArg{RooFit::ExternalConstraints(constraintPdfs)} {}
};

/// \see RooFit::FillColor()
struct FillColor : RooCmdArg {
   FillColor() = default;
   FillColor(TColorNumber color) : RooCmdArg{RooFit::FillColor(color)} {}
};

/// \see RooFit::FillStyle()
struct FillStyle : RooCmdArg {
   FillStyle() = default;
   FillStyle(Style_t style) : RooCmdArg{RooFit::FillStyle(style)} {}
   FillStyle(std::string const &style) : RooCmdArg{RooFit::FillStyle(style)} {}
};

/// \see RooFit::FitGauss()
struct FitGauss : RooCmdArg {
   FitGauss() = default;
   FitGauss(bool flag) : RooCmdArg{RooFit::FitGauss(flag)} {}
};

/// \see RooFit::FitModel()
struct FitModel : RooCmdArg {
   FitModel() = default;
   FitModel(RooAbsPdf &pdf) : RooCmdArg{RooFit::FitModel(pdf)} {}
};

/// \see RooFit::FixedPrecision()
struct FixedPrecision : RooCmdArg {
   FixedPrecision() = default;
   FixedPrecision(Int_t ndigit) : RooCmdArg{RooFit::FixedPrecision(ndigit)} {}
};

/// \see RooFit::FrameBins()
struct FrameBins : RooCmdArg {
   FrameBins() = default;
   FrameBins(Int_t nbins) : RooCmdArg{RooFit::FrameBins(nbins)} {}
};

/// \see RooFit::FrameRange()
struct FrameRange : RooCmdArg {
   FrameRange() = default;
   FrameRange(double xlo, double xhi) : RooCmdArg{RooFit::FrameRange(xlo, xhi)} {}
};

/// \see RooFit::GenBinned()
struct GenBinned : RooCmdArg {
   GenBinned() = default;
   GenBinned(const char *tag) : RooCmdArg{RooFit::GenBinned(tag)} {}
};

/// \see RooFit::GlobalObservables()
struct GlobalObservables : RooCmdArg {
   GlobalObservables() = default;
   GlobalObservables(const RooArgSet &args) : RooCmdArg{RooFit::GlobalObservables(args)} {}
};

/// \see RooFit::GlobalObservablesSource()
struct GlobalObservablesSource : RooCmdArg {
   GlobalObservablesSource() = default;
   GlobalObservablesSource(const char *sourceName) : RooCmdArg{RooFit::GlobalObservablesSource(sourceName)} {}
};

/// \see RooFit::GlobalObservablesTag()
struct GlobalObservablesTag : RooCmdArg {
   GlobalObservablesTag() = default;
   GlobalObservablesTag(const char *tagName) : RooCmdArg{RooFit::GlobalObservablesTag(tagName)} {}
};

/// \see RooFit::Hesse()
struct Hesse : RooCmdArg {
   Hesse() = default;
   Hesse(bool flag) : RooCmdArg{RooFit::Hesse(flag)} {}
};

/// \see RooFit::Import()
struct Import : RooCmdArg {
   Import() = default;
   Import(const char *state, TH1 &histo) : RooCmdArg{RooFit::Import(state, histo)} {}
   Import(const std::map<std::string, TH1 *> &arg) : RooCmdArg{RooFit::Import(arg)} {}
   Import(const char *state, RooDataHist &dhist) : RooCmdArg{RooFit::Import(state, dhist)} {}
   Import(const std::map<std::string, RooDataHist *> &arg) : RooCmdArg{RooFit::Import(arg)} {}
   Import(TH1 &histo, bool importDensity = false) : RooCmdArg{RooFit::Import(histo, importDensity)} {}
   Import(const char *state, RooAbsData &data) : RooCmdArg{RooFit::Import(state, data)} {}
   Import(const std::map<std::string, RooDataSet *> &arg) : RooCmdArg{RooFit::Import(arg)} {}
   Import(RooAbsData &data) : RooCmdArg{RooFit::Import(data)} {}
   Import(TTree &tree) : RooCmdArg{RooFit::Import(tree)} {}
};

/// \see RooFit::ImportFromFile()
struct ImportFromFile : RooCmdArg {
   ImportFromFile() = default;
   ImportFromFile(const char *fname, const char *tname) : RooCmdArg{RooFit::ImportFromFile(fname, tname)} {}
};

/// \see RooFit::Index()
struct Index : RooCmdArg {
   Index() = default;
   Index(RooCategory &icat) : RooCmdArg{RooFit::Index(icat)} {}
};

/// \see RooFit::InitialHesse()
struct InitialHesse : RooCmdArg {
   InitialHesse() = default;
   InitialHesse(bool flag) : RooCmdArg{RooFit::InitialHesse(flag)} {}
};

/// \see RooFit::IntegrateBins()
struct IntegrateBins : RooCmdArg {
   IntegrateBins() = default;
   IntegrateBins(double precision) : RooCmdArg{RooFit::IntegrateBins(precision)} {}
};

/// \see RooFit::IntrinsicBinning()
struct IntrinsicBinning : RooCmdArg {
   IntrinsicBinning() = default;
   IntrinsicBinning(bool flag) : RooCmdArg{RooFit::IntrinsicBinning(flag)} {}
};

/// \see RooFit::Invisible()
struct Invisible : RooCmdArg {
   Invisible() = default;
   Invisible(bool inv) : RooCmdArg{RooFit::Invisible(inv)} {}
};

/// \see RooFit::Label()
struct Label : RooCmdArg {
   Label() = default;
   Label(const char *str) : RooCmdArg{RooFit::Label(str)} {}
};

/// \see RooFit::LatexStyle()
struct LatexStyle : RooCmdArg {
   LatexStyle() = default;
   LatexStyle(bool flag) : RooCmdArg{RooFit::LatexStyle(flag)} {}
};

/// \see RooFit::LatexTableStyle()
struct LatexTableStyle : RooCmdArg {
   LatexTableStyle() = default;
   LatexTableStyle(bool flag) : RooCmdArg{RooFit::LatexTableStyle(flag)} {}
};

/// \see RooFit::Layout()
struct Layout : RooCmdArg {
   Layout() = default;
   Layout(double xmin, double xmax = 0.99, double ymin = 0.95) : RooCmdArg{RooFit::Layout(xmin, xmax, ymin)} {}
};

/// \see RooFit::LineColor()
struct LineColor : RooCmdArg {
   LineColor() = default;
   LineColor(TColorNumber color) : RooCmdArg{RooFit::LineColor(color)} {}
};

/// \see RooFit::LineStyle()
struct LineStyle : RooCmdArg {
   LineStyle() = default;
   LineStyle(Style_t style) : RooCmdArg{RooFit::LineStyle(style)} {}
   LineStyle(std::string const &style) : RooCmdArg{RooFit::LineStyle(style)} {}
};

/// \see RooFit::LineWidth()
struct LineWidth : RooCmdArg {
   LineWidth() = default;
   LineWidth(Width_t width) : RooCmdArg{RooFit::LineWidth(width)} {}
};

/// \see RooFit::Link()
struct Link : RooCmdArg {
   Link() = default;
   Link(const char *state, RooAbsData &data) : RooCmdArg{RooFit::Link(state, data)} {}
   Link(const std::map<std::string, RooAbsData *> &arg) : RooCmdArg{RooFit::Link(arg)} {}
};

/// \see RooFit::MarkerColor()
struct MarkerColor : RooCmdArg {
   MarkerColor() = default;
   MarkerColor(TColorNumber color) : RooCmdArg{RooFit::MarkerColor(color)} {}
};

/// \see RooFit::MarkerSize()
struct MarkerSize : RooCmdArg {
   MarkerSize() = default;
   MarkerSize(Size_t size) : RooCmdArg{RooFit::MarkerSize(size)} {}
};

/// \see RooFit::MarkerStyle()
struct MarkerStyle : RooCmdArg {
   MarkerStyle() = default;
   MarkerStyle(Style_t style) : RooCmdArg{RooFit::MarkerStyle(style)} {}
   MarkerStyle(std::string const &style) : RooCmdArg{RooFit::MarkerStyle(style)} {}
};

/// \see RooFit::MaxCalls()
struct MaxCalls : RooCmdArg {
   MaxCalls() = default;
   MaxCalls(int n) : RooCmdArg{RooFit::MaxCalls(n)} {}
};

/// \see RooFit::Minimizer()
struct Minimizer : RooCmdArg {
   Minimizer() = default;
   Minimizer(const char *type, const char *alg = nullptr) : RooCmdArg{RooFit::Minimizer(type, alg)} {}
};

/// \see RooFit::Minos()
struct Minos : RooCmdArg {
   Minos() = default;
   Minos(bool flag) : RooCmdArg{RooFit::Minos(flag)} {}
   Minos(const RooArgSet &minosArgs) : RooCmdArg{RooFit::Minos(minosArgs)} {}
};

/// \see RooFit::ModularL()
struct ModularL : RooCmdArg {
   ModularL() = default;
   ModularL(bool flag) : RooCmdArg{RooFit::ModularL(flag)} {}
};

/// \see RooFit::MoveToBack()
struct MoveToBack : RooCmdArg {
   MoveToBack() = default;
   MoveToBack(bool flag) : RooCmdArg{flag ? RooFit::MoveToBack() : RooCmdArg{}} {}
};

/// \see RooFit::Name()
struct Name : RooCmdArg {
   Name() = default;
   Name(const char *name) : RooCmdArg{RooFit::Name(name)} {}
};

/// \see RooFit::NoRecursion()
struct NoRecursion : RooCmdArg {
   NoRecursion() = default;
   NoRecursion(bool flag) : RooCmdArg{RooFit::NoRecursion(flag)} {}
};

/// \see RooFit::NormRange()
struct NormRange : RooCmdArg {
   NormRange() = default;
   NormRange(const char *rangeNameList) : RooCmdArg{RooFit::NormRange(rangeNameList)} {}
};

/// \see RooFit::NormSet()
struct NormSet : RooCmdArg {
   NormSet() = default;
   NormSet(const RooArgSet &args) : RooCmdArg{RooFit::NormSet(args)} {}
};

/// \see RooFit::Normalization()
struct Normalization : RooCmdArg {
   Normalization() = default;
   Normalization(double scaleFactor) : RooCmdArg{RooFit::Normalization(scaleFactor)} {}
   Normalization(double scaleFactor, Int_t scaleType) : RooCmdArg{RooFit::Normalization(scaleFactor, scaleType)} {}
};

/// \see RooFit::NumCPU()
struct NumCPU : RooCmdArg {
   NumCPU() = default;
   NumCPU(Int_t nCPU, Int_t interleave = 0) : RooCmdArg{RooFit::NumCPU(nCPU, interleave)} {}
};

/// \see RooFit::NumEvents()
struct NumEvents : RooCmdArg {
   NumEvents() = default;
   NumEvents(Int_t numEvents) : RooCmdArg{RooFit::NumEvents(numEvents)} {}
   NumEvents(double numEvents) : RooCmdArg{RooFit::NumEvents(numEvents)} {}
};

/// \see RooFit::NumIntConfig()
struct NumIntConfig : RooCmdArg {
   NumIntConfig() = default;
   NumIntConfig(const RooNumIntConfig &cfg) : RooCmdArg{RooFit::NumIntConfig(cfg)} {}
};

/// \see RooFit::ObjectName()
struct ObjectName : RooCmdArg {
   ObjectName() = default;
   ObjectName(const char *name) : RooCmdArg{RooFit::ObjectName(name)} {}
};

/// \see RooFit::Offset()
struct Offset : RooCmdArg {
   Offset() = default;
   Offset(std::string const &mode) : RooCmdArg{RooFit::Offset(mode)} {}
   Offset(const char *mode) : RooCmdArg{RooFit::Offset(mode)} {}
   Offset(bool flag) : RooCmdArg{RooFit::Offset(flag)} {}
};

/// \see RooFit::OutputFile()
struct OutputFile : RooCmdArg {
   OutputFile() = default;
   OutputFile(const char *fileName) : RooCmdArg{RooFit::OutputFile(fileName)} {}
};

/// \see RooFit::OutputStream()
struct OutputStream : RooCmdArg {
   OutputStream() = default;
   OutputStream(std::ostream &os) : RooCmdArg{RooFit::OutputStream(os)} {}
};

/// \see RooFit::OwnLinked()
struct OwnLinked : RooCmdArg {
   OwnLinked() = default;
   OwnLinked(bool flag) : RooCmdArg{flag ? RooFit::OwnLinked() : RooCmdArg{}} {}
};

/// \see RooFit::Experimental::ParallelDescentOptions()
struct ParallelDescentOptions : RooCmdArg {
   ParallelDescentOptions() = default;
   ParallelDescentOptions(bool enable, int splitStrategy = 0, int numSplits = 0)
      : RooCmdArg{RooFit::Experimental::ParallelDescentOptions(enable, splitStrategy, numSplits)}
   {
   }
};

/// \see RooFit::Experimental::ParallelGradientOptions()
struct ParallelGradientOptions : RooCmdArg {
   ParallelGradientOptions() = default;
   ParallelGradientOptions(bool enable, int orderStrategy = 0, int chainFactor = 1)
      : RooCmdArg{RooFit::Experimental::ParallelGradientOptions(enable, orderStrategy, chainFactor)}
   {
   }
};

/// \see RooFit::Parallelize()
struct Parallelize : RooCmdArg {
   Parallelize() = default;
   Parallelize(int nWorkers) : RooCmdArg{RooFit::Parallelize(nWorkers)} {}
};

/// \see RooFit::Parameters()
struct Parameters : RooCmdArg {
   Parameters() = default;
   Parameters(const RooArgSet &params) : RooCmdArg{RooFit::Parameters(params)} {}
};

/// \see RooFit::Precision()
struct Precision : RooCmdArg {
   Precision() = default;
   Precision(double prec) : RooCmdArg{RooFit::Precision(prec)} {}
};

/// \see RooFit::PrefitDataFraction()
struct PrefitDataFraction : RooCmdArg {
   PrefitDataFraction() = default;
   PrefitDataFraction(double data_ratio) : RooCmdArg{RooFit::PrefitDataFraction(data_ratio)} {}
};

/// \see RooFit::Prefix()
struct Prefix : RooCmdArg {
   Prefix() = default;
   Prefix(bool flag) : RooCmdArg{RooFit::Prefix(flag)} {}
};

/// \see RooFit::PrintEvalErrors()
struct PrintEvalErrors : RooCmdArg {
   PrintEvalErrors() = default;
   PrintEvalErrors(Int_t numErrors) : RooCmdArg{RooFit::PrintEvalErrors(numErrors)} {}
};

/// \see RooFit::PrintLevel()
struct PrintLevel : RooCmdArg {
   PrintLevel() = default;
   PrintLevel(Int_t code) : RooCmdArg{RooFit::PrintLevel(code)} {}
};

/// \see RooFit::ProjWData()
struct ProjWData : RooCmdArg {
   ProjWData() = default;
   ProjWData(const RooAbsData &projData, bool binData = false) : RooCmdArg{RooFit::ProjWData(projData, binData)} {}
   ProjWData(const RooArgSet &projSet, const RooAbsData &projData, bool binData = false)
      : RooCmdArg{RooFit::ProjWData(projSet, projData, binData)}
   {
   }
};

/// \see RooFit::Project()
struct Project : RooCmdArg {
   Project() = default;
   Project(const RooArgSet &projSet) : RooCmdArg{RooFit::Project(projSet)} {}
};

/// \see RooFit::ProjectedObservables()
struct ProjectedObservables : RooCmdArg {
   ProjectedObservables() = default;
   ProjectedObservables(const RooArgSet &args) : RooCmdArg{RooFit::ProjectedObservables(args)} {}
};

/// \see RooFit::ProjectionRange()
struct ProjectionRange : RooCmdArg {
   ProjectionRange() = default;
   ProjectionRange(const char *rangeName) : RooCmdArg{RooFit::ProjectionRange(rangeName)} {}
};

/// \see RooFit::ProtoData()
struct ProtoData : RooCmdArg {
   ProtoData() = default;
   ProtoData(const RooAbsData &protoData, bool randomizeOrder = false, bool resample = false)
      : RooCmdArg{RooFit::ProtoData(protoData, randomizeOrder, resample)}
   {
   }
};

/// \see RooFit::Range()
struct Range : RooCmdArg {
   Range() = default;
   Range(const char *rangeName, bool adjustNorm = true) : RooCmdArg{RooFit::Range(rangeName, adjustNorm)} {}
   Range(double lo, double hi, bool adjustNorm = true) : RooCmdArg{RooFit::Range(lo, hi, adjustNorm)} {}
};

/// \see RooFit::RecoverFromUndefinedRegions()
struct RecoverFromUndefinedRegions : RooCmdArg {
   RecoverFromUndefinedRegions() = default;
   RecoverFromUndefinedRegions(double strength) : RooCmdArg{RooFit::RecoverFromUndefinedRegions(strength)} {}
};

/// \see RooFit::RecycleConflictNodes()
struct RecycleConflictNodes : RooCmdArg {
   RecycleConflictNodes() = default;
   RecycleConflictNodes(bool flag) : RooCmdArg{RooFit::RecycleConflictNodes(flag)} {}
};

/// \see RooFit::RefreshNorm()
struct RefreshNorm : RooCmdArg {
   RefreshNorm() = default;
   RefreshNorm(bool flag) : RooCmdArg{flag ? RooFit::RefreshNorm() : RooCmdArg{}} {}
};

/// \see RooFit::Rename()
struct Rename : RooCmdArg {
   Rename() = default;
   Rename(const char *suffix) : RooCmdArg{RooFit::Rename(suffix)} {}
};

/// \see RooFit::RenameAllNodes()
struct RenameAllNodes : RooCmdArg {
   RenameAllNodes() = default;
   RenameAllNodes(const char *suffix) : RooCmdArg{RooFit::RenameAllNodes(suffix)} {}
};

/// \see RooFit::RenameAllVariables()
struct RenameAllVariables : RooCmdArg {
   RenameAllVariables() = default;
   RenameAllVariables(const char *suffix) : RooCmdArg{RooFit::RenameAllVariables(suffix)} {}
};

/// \see RooFit::RenameAllVariablesExcept()
struct RenameAllVariablesExcept : RooCmdArg {
   RenameAllVariablesExcept() = default;
   RenameAllVariablesExcept(const char *suffix, const char *exceptionList)
      : RooCmdArg{RooFit::RenameAllVariablesExcept(suffix, exceptionList)}
   {
   }
};

/// \see RooFit::RenameConflictNodes()
struct RenameConflictNodes : RooCmdArg {
   RenameConflictNodes() = default;
   RenameConflictNodes(const char *suffix, bool renameOrigNodes = false)
      : RooCmdArg{RooFit::RenameConflictNodes(suffix, renameOrigNodes)}
   {
   }
};

/// \see RooFit::RenameVariable()
struct RenameVariable : RooCmdArg {
   RenameVariable() = default;
   RenameVariable(const char *inputName, const char *outputName)
      : RooCmdArg{RooFit::RenameVariable(inputName, outputName)}
   {
   }
};

/// \see RooFit::Rescale()
struct Rescale : RooCmdArg {
   Rescale() = default;
   Rescale(double factor) : RooCmdArg{RooFit::Rescale(factor)} {}
};

/// \see RooFit::Restrict()
struct Restrict : RooCmdArg {
   Restrict() = default;
   Restrict(const char *catName, const char *stateNameList) : RooCmdArg{RooFit::Restrict(catName, stateNameList)} {}
};

/// \see RooFit::Save()
struct Save : RooCmdArg {
   Save() = default;
   Save(bool flag) : RooCmdArg{RooFit::Save(flag)} {}
};

/// \see RooFit::Scaling()
struct Scaling : RooCmdArg {
   Scaling() = default;
   Scaling(bool flag) : RooCmdArg{RooFit::Scaling(flag)} {}
};

/// \see RooFit::ScanAllCdf()
struct ScanAllCdf : RooCmdArg {
   ScanAllCdf() = default;
   ScanAllCdf(bool flag) : RooCmdArg{flag ? RooFit::ScanAllCdf() : RooCmdArg{}} {}
};

/// \see RooFit::ScanNoCdf()
struct ScanNoCdf : RooCmdArg {
   ScanNoCdf() = default;
   ScanNoCdf(bool flag) : RooCmdArg{flag ? RooFit::ScanNoCdf() : RooCmdArg{}} {}
};

/// \see RooFit::ScanNumCdf()
struct ScanNumCdf : RooCmdArg {
   ScanNumCdf() = default;
   ScanNumCdf(bool flag) : RooCmdArg{flag ? RooFit::ScanNumCdf() : RooCmdArg{}} {}
};

/// \see RooFit::ScanParameters()
struct ScanParameters : RooCmdArg {
   ScanParameters() = default;
   ScanParameters(Int_t nbins, Int_t intOrder) : RooCmdArg{RooFit::ScanParameters(nbins, intOrder)} {}
};

/// \see RooFit::SelectVars()
struct SelectVars : RooCmdArg {
   SelectVars() = default;
   SelectVars(const RooArgSet &vars) : RooCmdArg{RooFit::SelectVars(vars)} {}
};

/// \see RooFit::ShiftToZero()
struct ShiftToZero : RooCmdArg {
   ShiftToZero() = default;
   ShiftToZero(bool flag) : RooCmdArg{flag ? RooFit::ShiftToZero() : RooCmdArg{}} {}
};

/// \see RooFit::ShowConstants()
struct ShowConstants : RooCmdArg {
   ShowConstants() = default;
   ShowConstants(bool flag) : RooCmdArg{RooFit::ShowConstants(flag)} {}
};

/// \see RooFit::ShowProgress()
struct ShowProgress : RooCmdArg {
   ShowProgress() = default;
   ShowProgress(bool flag) : RooCmdArg{flag ? RooFit::ShowProgress() : RooCmdArg{}} {}
};

/// \see RooFit::Sibling()
struct Sibling : RooCmdArg {
   Sibling() = default;
   Sibling(const RooAbsCollection &sibling) : RooCmdArg{RooFit::Sibling(sibling)} {}
};

/// \see RooFit::Silence()
struct Silence : RooCmdArg {
   Silence() = default;
   Silence(bool flag) : RooCmdArg{RooFit::Silence(flag)} {}
};

/// \see RooFit::Slice()
struct Slice : RooCmdArg {
   Slice() = default;
   Slice(const RooArgSet &sliceSet) : RooCmdArg{RooFit::Slice(sliceSet)} {}
   Slice(RooCategory &cat, const char *label) : RooCmdArg{RooFit::Slice(cat, label)} {}
   Slice(std::map<RooCategory *, std::string> const &arg) : RooCmdArg{RooFit::Slice(arg)} {}
};

/// \see RooFit::SplitParam()
struct SplitParam : RooCmdArg {
   SplitParam() = default;
   SplitParam(const char *varname, const char *catname) : RooCmdArg{RooFit::SplitParam(varname, catname)} {}
   SplitParam(const RooRealVar &var, const RooAbsCategory &cat) : RooCmdArg{RooFit::SplitParam(var, cat)} {}
};

/// \see RooFit::SplitParamConstrained()
struct SplitParamConstrained : RooCmdArg {
   SplitParamConstrained() = default;
   SplitParamConstrained(const char *varname, const char *catname, const char *rsname)
      : RooCmdArg{RooFit::SplitParamConstrained(varname, catname, rsname)}
   {
   }
   SplitParamConstrained(const RooRealVar &var, const RooAbsCategory &cat, const char *rsname)
      : RooCmdArg{RooFit::SplitParamConstrained(var, cat, rsname)}
   {
   }
};

/// \see RooFit::SplitRange()
struct SplitRange : RooCmdArg {
   SplitRange() = default;
   SplitRange(bool flag) : RooCmdArg{RooFit::SplitRange(flag)} {}
};

/// \see RooFit::StoreAsymError()
struct StoreAsymError : RooCmdArg {
   StoreAsymError() = default;
   StoreAsymError(const RooArgSet &aset) : RooCmdArg{RooFit::StoreAsymError(aset)} {}
};

/// \see RooFit::StoreError()
struct StoreError : RooCmdArg {
   StoreError() = default;
   StoreError(const RooArgSet &aset) : RooCmdArg{RooFit::StoreError(aset)} {}
};

/// \see RooFit::Strategy()
struct Strategy : RooCmdArg {
   Strategy() = default;
   Strategy(Int_t code) : RooCmdArg{RooFit::Strategy(code)} {}
};

/// \see RooFit::SumCoefRange()
struct SumCoefRange : RooCmdArg {
   SumCoefRange() = default;
   SumCoefRange(const char *rangeName) : RooCmdArg{RooFit::SumCoefRange(rangeName)} {}
};

/// \see RooFit::SumW2Error()
struct SumW2Error : RooCmdArg {
   SumW2Error() = default;
   SumW2Error(bool flag) : RooCmdArg{RooFit::SumW2Error(flag)} {}
};

/// \see RooFit::SupNormSet()
struct SupNormSet : RooCmdArg {
   SupNormSet() = default;
   SupNormSet(const RooArgSet &nset) : RooCmdArg{RooFit::SupNormSet(nset)} {}
};

/// \see RooFit::TLatexStyle()
struct TLatexStyle : RooCmdArg {
   TLatexStyle() = default;
   TLatexStyle(bool flag) : RooCmdArg{RooFit::TLatexStyle(flag)} {}
};

/// \see RooFit::TagName()
struct TagName : RooCmdArg {
   TagName() = default;
   TagName(const char *name) : RooCmdArg{RooFit::TagName(name)} {}
};

/// \see RooFit::Timer()
struct Timer : RooCmdArg {
   Timer() = default;
   Timer(bool flag) : RooCmdArg{RooFit::Timer(flag)} {}
};

/// \see RooFit::TimingAnalysis()
struct TimingAnalysis : RooCmdArg {
   TimingAnalysis() = default;
   TimingAnalysis(bool timingAnalysis) : RooCmdArg{RooFit::TimingAnalysis(timingAnalysis)} {}
};

/// \see RooFit::Title()
struct Title : RooCmdArg {
   Title() = default;
   Title(const char *name) : RooCmdArg{RooFit::Title(name)} {}
};

/// \see RooFit::Topic()
struct Topic : RooCmdArg {
   Topic() = default;
   Topic(Int_t topic) : RooCmdArg{RooFit::Topic(topic)} {}
};

/// \see RooFit::VLines()
struct VLines : RooCmdArg {
   VLines() = default;
   VLines(bool flag) : RooCmdArg{flag ? RooFit::VLines() : RooCmdArg{}} {}
};

/// \see RooFit::VerbatimName()
struct VerbatimName : RooCmdArg {
   VerbatimName() = default;
   VerbatimName(bool flag) : RooCmdArg{RooFit::VerbatimName(flag)} {}
};

/// \see RooFit::Verbose()
struct Verbose : RooCmdArg {
   Verbose() = default;
   Verbose(bool flag) : RooCmdArg{RooFit::Verbose(flag)} {}
};

/// \see RooFit::VisualizeError()
struct VisualizeError : RooCmdArg {
   VisualizeError() = default;
   VisualizeError(const RooDataSet &paramData, double Z = 1) : RooCmdArg{RooFit::VisualizeError(paramData, Z)} {}
   VisualizeError(const RooFitResult &fitres, double Z = 1, bool linearMethod = true)
      : RooCmdArg{RooFit::VisualizeError(fitres, Z, linearMethod)}
   {
   }
   VisualizeError(const RooFitResult &fitres, const RooArgSet &param, double Z = 1, bool linearMethod = true)
      : RooCmdArg{RooFit::VisualizeError(fitres, param, Z, linearMethod)}
   {
   }
};

/// \see RooFit::Warnings()
struct Warnings : RooCmdArg {
   Warnings() = default;
   Warnings(bool flag) : RooCmdArg{RooFit::Warnings(flag)} {}
};

/// \see RooFit::Weight()
struct Weight : RooCmdArg {
   Weight() = default;
   Weight(double wgt) : RooCmdArg{RooFit::Weight(wgt)} {}
};

/// \see RooFit::WeightVar()
struct WeightVar : RooCmdArg {
   WeightVar() = default;
   WeightVar(const char *name, bool reinterpretAsWeight = false)
      : RooCmdArg{RooFit::WeightVar(name, reinterpretAsWeight)}
   {
   }
   WeightVar(const RooRealVar &arg, bool reinterpretAsWeight = false)
      : RooCmdArg{RooFit::WeightVar(arg, reinterpretAsWeight)}
   {
   }
};

/// \see RooFit::What()
struct What : RooCmdArg {
   What() = default;
   What(const char *str) : RooCmdArg{RooFit::What(str)} {}
};

/// \see RooFit::XErrorSize()
struct XErrorSize : RooCmdArg {
   XErrorSize() = default;
   XErrorSize(double width) : RooCmdArg{RooFit::XErrorSize(width)} {}
};

} // namespace Kwarg
} // namespace Detail

/// Keyword arguments for RooAbsPdf::fitTo().
/// Each member corresponds to the RooFit command argument of the same name.
struct FitToKwargs {
   RooFit::Detail::Kwarg::AsymptoticError AsymptoticError;
   RooFit::Detail::Kwarg::CloneData CloneData;
   RooFit::Detail::Kwarg::ConditionalObservables ConditionalObservables;
   RooFit::Detail::Kwarg::Constrain Constrain;
   RooFit::Detail::Kwarg::EvalBackend EvalBackend;
   RooFit::Detail::Kwarg::EvalErrorWall EvalErrorWall;
   RooFit::Detail::Kwarg::Extended Extended;
   RooFit::Detail::Kwarg::ExternalConstraints ExternalConstraints;
   RooFit::Detail::Kwarg::GlobalObservables GlobalObservables;
   RooFit::Detail::Kwarg::GlobalObservablesSource GlobalObservablesSource;
   RooFit::Detail::Kwarg::GlobalObservablesTag GlobalObservablesTag;
   RooFit::Detail::Kwarg::Hesse Hesse;
   RooFit::Detail::Kwarg::InitialHesse InitialHesse;
   RooFit::Detail::Kwarg::IntegrateBins IntegrateBins;
   RooFit::Detail::Kwarg::MaxCalls MaxCalls;
   RooFit::Detail::Kwarg::Minimizer Minimizer;
   RooFit::Detail::Kwarg::Minos Minos;
   RooFit::Detail::Kwarg::ModularL ModularL;
   RooFit::Detail::Kwarg::NumCPU NumCPU;
   RooFit::Detail::Kwarg::Offset Offset;
   RooFit::Detail::Kwarg::ParallelDescentOptions ParallelDescentOptions;
   RooFit::Detail::Kwarg::ParallelGradientOptions ParallelGradientOptions;
   RooFit::Detail::Kwarg::Parallelize Parallelize;
   RooFit::Detail::Kwarg::PrefitDataFraction PrefitDataFraction;
   RooFit::Detail::Kwarg::PrintEvalErrors PrintEvalErrors;
   RooFit::Detail::Kwarg::PrintLevel PrintLevel;
   RooFit::Detail::Kwarg::ProjectedObservables ProjectedObservables;
   RooFit::Detail::Kwarg::Range Range;
   RooFit::Detail::Kwarg::RecoverFromUndefinedRegions RecoverFromUndefinedRegions;
   RooFit::Detail::Kwarg::Save Save;
   RooFit::Detail::Kwarg::SplitRange SplitRange;
   RooFit::Detail::Kwarg::Strategy Strategy;
   RooFit::Detail::Kwarg::SumCoefRange SumCoefRange;
   RooFit::Detail::Kwarg::SumW2Error SumW2Error;
   RooFit::Detail::Kwarg::Timer Timer;
   RooFit::Detail::Kwarg::TimingAnalysis TimingAnalysis;
   RooFit::Detail::Kwarg::Verbose Verbose;
   RooFit::Detail::Kwarg::Warnings Warnings;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&AsymptoticError, &CloneData, &ConditionalObservables, &Constrain, &EvalBackend, &EvalErrorWall,
              &Extended, &ExternalConstraints, &GlobalObservables, &GlobalObservablesSource, &GlobalObservablesTag,
              &Hesse, &InitialHesse, &IntegrateBins, &MaxCalls, &Minimizer, &Minos, &ModularL, &NumCPU, &Offset,
              &ParallelDescentOptions, &ParallelGradientOptions, &Parallelize, &PrefitDataFraction, &PrintEvalErrors,
              &PrintLevel, &ProjectedObservables, &Range, &RecoverFromUndefinedRegions, &Save, &SplitRange, &Strategy,
              &SumCoefRange, &SumW2Error, &Timer, &TimingAnalysis, &Verbose, &Warnings};
      // clang-format on
   }
};

/// Keyword arguments for RooFit::Format().
/// Each member corresponds to the RooFit command argument of the same name.
struct FormatKwargs {
   RooFit::Detail::Kwarg::AutoPrecision AutoPrecision;
   RooFit::Detail::Kwarg::FixedPrecision FixedPrecision;
   RooFit::Detail::Kwarg::LatexStyle LatexStyle;
   RooFit::Detail::Kwarg::LatexTableStyle LatexTableStyle;
   RooFit::Detail::Kwarg::TLatexStyle TLatexStyle;
   RooFit::Detail::Kwarg::VerbatimName VerbatimName;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&AutoPrecision, &FixedPrecision, &LatexStyle, &LatexTableStyle, &TLatexStyle, &VerbatimName};
      // clang-format on
   }
};

/// Keyword arguments for RooAbsRealLValue::frame().
/// Each member corresponds to the RooFit command argument of the same name.
struct FrameKwargs {
   RooFit::Detail::Kwarg::AutoRange AutoRange;
   RooFit::Detail::Kwarg::AutoSymRange AutoSymRange;
   RooFit::Detail::Kwarg::Bins Bins;
   RooFit::Detail::Kwarg::Name Name;
   RooFit::Detail::Kwarg::Range Range;
   RooFit::Detail::Kwarg::Title Title;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&AutoRange, &AutoSymRange, &Bins, &Name, &Range, &Title};
      // clang-format on
   }
};

/// Keyword arguments for RooFit::YVar().
/// Each member corresponds to the RooFit command argument of the same name.
struct YVarKwargs {
   RooFit::Detail::Kwarg::AutoBinning AutoBinning;
   RooFit::Detail::Kwarg::AutoSymBinning AutoSymBinning;
   RooFit::Detail::Kwarg::Binning Binning;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&AutoBinning, &AutoSymBinning, &Binning};
      // clang-format on
   }
};

/// Keyword arguments for RooFit::ZVar().
/// Each member corresponds to the RooFit command argument of the same name.
struct ZVarKwargs {
   RooFit::Detail::Kwarg::AutoBinning AutoBinning;
   RooFit::Detail::Kwarg::AutoSymBinning AutoSymBinning;
   RooFit::Detail::Kwarg::Binning Binning;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&AutoBinning, &AutoSymBinning, &Binning};
      // clang-format on
   }
};

RooCmdArg FitOptions(FitToKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                     const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                     const RooCmdArg &arg6 = {});
RooCmdArg Format(const char *what, FormatKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                 const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                 const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {});
RooCmdArg Frame(FrameKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
                const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
                const RooCmdArg &arg6 = {});
RooCmdArg YVar(const RooAbsRealLValue &var, YVarKwargs const &kwargs, const RooCmdArg &arg = {});
RooCmdArg ZVar(const RooAbsRealLValue &var, ZVarKwargs const &kwargs, const RooCmdArg &arg = {});

namespace Detail {
namespace Kwarg {

/// \see RooFit::FitOptions()
struct FitOptions : RooCmdArg {
   FitOptions() = default;
   FitOptions(const RooCmdArg &arg1, const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
              const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {})
      : RooCmdArg{RooFit::FitOptions(arg1, arg2, arg3, arg4, arg5, arg6)}
   {
   }
   FitOptions(RooFit::FitToKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
              const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
              const RooCmdArg &arg6 = {})
      : RooCmdArg{RooFit::FitOptions(kwargs, arg1, arg2, arg3, arg4, arg5, arg6)}
   {
   }
};

/// \see RooFit::Format()
struct Format : RooCmdArg {
   Format() = default;
   Format(const char *what, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {},
          const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {},
          const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
      : RooCmdArg{RooFit::Format(what, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)}
   {
   }
   Format(const char *what, RooFit::FormatKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
          const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {},
          const RooCmdArg &arg6 = {}, const RooCmdArg &arg7 = {}, const RooCmdArg &arg8 = {})
      : RooCmdArg{RooFit::Format(what, kwargs, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)}
   {
   }
};

/// \see RooFit::Frame()
struct Frame : RooCmdArg {
   Frame() = default;
   Frame(const RooCmdArg &arg1, const RooCmdArg &arg2 = {}, const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {},
         const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {})
      : RooCmdArg{RooFit::Frame(arg1, arg2, arg3, arg4, arg5, arg6)}
   {
   }
   Frame(RooFit::FrameKwargs const &kwargs, const RooCmdArg &arg1 = {}, const RooCmdArg &arg2 = {},
         const RooCmdArg &arg3 = {}, const RooCmdArg &arg4 = {}, const RooCmdArg &arg5 = {}, const RooCmdArg &arg6 = {})
      : RooCmdArg{RooFit::Frame(kwargs, arg1, arg2, arg3, arg4, arg5, arg6)}
   {
   }
};

/// \see RooFit::YVar()
struct YVar : RooCmdArg {
   YVar() = default;
   YVar(const RooAbsRealLValue &var, const RooCmdArg &arg = {}) : RooCmdArg{RooFit::YVar(var, arg)} {}
   YVar(const RooAbsRealLValue &var, RooFit::YVarKwargs const &kwargs, const RooCmdArg &arg = {})
      : RooCmdArg{RooFit::YVar(var, kwargs, arg)}
   {
   }
};

/// \see RooFit::ZVar()
struct ZVar : RooCmdArg {
   ZVar() = default;
   ZVar(const RooAbsRealLValue &var, const RooCmdArg &arg = {}) : RooCmdArg{RooFit::ZVar(var, arg)} {}
   ZVar(const RooAbsRealLValue &var, RooFit::ZVarKwargs const &kwargs, const RooCmdArg &arg = {})
      : RooCmdArg{RooFit::ZVar(var, kwargs, arg)}
   {
   }
};

} // namespace Kwarg
} // namespace Detail
} // namespace RooFit

/// Keyword arguments for RooAbsCollection::printLatex().
/// Each member corresponds to the RooFit command argument of the same name.
struct RooAbsCollection::PrintLatexKwargs {
   RooFit::Detail::Kwarg::Columns Columns;
   RooFit::Detail::Kwarg::Format Format;
   RooFit::Detail::Kwarg::OutputFile OutputFile;
   RooFit::Detail::Kwarg::Sibling Sibling;

   std::vector<RooCmdArg const *> cmdArgs() const
   {
      // clang-format off
      return {&Columns, &Format, &OutputFile, &Sibling};
      // clang-format on
   }
};

#endif
