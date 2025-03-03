// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
// Author: Giovanni Petrucciani (UCSD) (log-interpolation)
// Modified by Hongtao Yang for xmlAnaWSBuilder
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//_________________________________________________
/*
BEGIN_HTML
<p>
</p>
END_HTML
*/
//

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooTrace.h"

#include "TMath.h"

#include "ResponseFunction.h"

ClassImp(ResponseFunction)

//_____________________________________________________________________________
ResponseFunction::ResponseFunction()
{
  // Default constructor
  _interpBoundary = 1.;
  _logInit = kFALSE;
  TRACE_CREATE
}

//_____________________________________________________________________________
ResponseFunction::ResponseFunction(const char *name, const char *title,
                                   const RooArgList &paramList,
                                   std::vector<double> nominal, const RooArgList &lowList, const RooArgList &highList, std::vector<int> code)
    : RooAbsReal(name, title),
      _nominal(nominal),
      _paramList("paramList", "List of paramficients", this),
      _lowList("lowList", "List of low variation parameters", this),
      _highList("highList", "List of high variation parameters", this),
      _interpCode(code),
      _interpBoundary(1.)
{

  _logInit = kFALSE;

  for (auto* it : paramList)
  {
    if (!dynamic_cast<RooAbsReal *>(it))
    {
      coutE(InputArguments) << "ResponseFunction::ctor(" << GetName() << ") ERROR: paramficient " << (it)->GetName()
                            << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _paramList.add(*it);
  }

  for (auto* it : lowList)
  {
    if (!dynamic_cast<RooAbsReal *>(it))
    {
      coutE(InputArguments) << "ResponseFunction::ctor(" << GetName() << ") ERROR: paramficient " << (it)->GetName()
                            << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _lowList.add(*it);
  }

  for (auto* it : highList )
  {
    if (!dynamic_cast<RooAbsReal *>(it))
    {
      coutE(InputArguments) << "ResponseFunction::ctor(" << GetName() << ") ERROR: paramficient " << (it)->GetName()
                            << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _highList.add(*it);
  }

  TRACE_CREATE
}

//_____________________________________________________________________________
ResponseFunction::ResponseFunction(const ResponseFunction &other, const char *name)
    : RooAbsReal(other, name),
      _nominal(other._nominal),
      _paramList("paramList", this, other._paramList),
      _lowList("lowList", this, other._lowList),
      _highList("highList", this, other._highList),
      _interpCode(other._interpCode), _interpBoundary(other._interpBoundary)

{
  // Copy constructor
  _logInit = kFALSE;
  TRACE_CREATE
}

//_____________________________________________________________________________
ResponseFunction::~ResponseFunction()
{
  // Destructor

  TRACE_DESTROY
}

//_____________________________________________________________________________
void ResponseFunction::setInterpCode(RooAbsReal &param, int code)
{

  int index = _paramList.index(&param);
  if (index < 0)
  {
    coutE(InputArguments) << "ResponseFunction::setInterpCode ERROR:  " << param.GetName()
                          << " is not in list" << std::endl;
  }
  else
  {
    coutW(InputArguments) << "ResponseFunction::setInterpCode :  " << param.GetName()
                          << " is now " << code << std::endl;
    _interpCode.at(index) = code;
  }
  // GHL: Adding suggestion by Swagato:
  _logInit = kFALSE;
  setValueDirty();
}

//_____________________________________________________________________________
void ResponseFunction::setAllInterpCodes(int code)
{

  for (unsigned int i = 0; i < _interpCode.size(); ++i)
  {
    _interpCode.at(i) = code;
  }
  // GHL: Adding suggestion by Swagato:
  _logInit = kFALSE;
  setValueDirty();
}

//_____________________________________________________________________________
double ResponseFunction::PolyInterpValue(int i, double x) const
{
  // If this is running for first time
  if (!_logInit)
  {
    // Create polynomial array
    const unsigned int n = _lowList.getSize();
    assert(n == (unsigned int)_highList.getSize());
    _polCoeff.resize(n * 6);
    // cache the polynomial coefficient values
    // which do not dpened on x but on the boundaries values
    for (unsigned j = 0; j < n; j++)
      cacheCoef(j);
    _logInit = kTRUE;
  }
  else
  {
    // Check upper and lower variations of i to see whether it is constant
    // If not, need to update the coefficient everytime
    if (!_highList[i].isConstant() || !_lowList[i].isConstant())
      cacheCoef(i);
  }

  // GHL: Swagato's suggestions
  // if( _low[i] == 0 ) _low[i] = 0.0001;
  // if( _high[i] == 0 ) _high[i] = 0.0001;

  // get pointer to location of coefficients in the vector
  const double *coeff = &_polCoeff.front() + 6 * i;

  double a = coeff[0];
  double b = coeff[1];
  double c = coeff[2];
  double d = coeff[3];
  double e = coeff[4];
  double f = coeff[5];

  // evaluate the 6-th degree polynomial using Horner's method
  double value = 1. + x * (a + x * (b + x * (c + x * (d + x * (e + x * f)))));
  return value;
}

//_____________________________________________________________________________
Double_t ResponseFunction::evaluate() const
{
  // Calculate and return value of polynomial

  Double_t total = 1;

  for (int i=0; i < _paramList.getSize(); i++) {
    const auto& param = static_cast<RooAbsReal&>(_paramList[i]);
    const auto& low = static_cast<RooAbsReal&>(_lowList[i]);
    const auto& high = static_cast<RooAbsReal&>(_highList[i]);

    Int_t icode = _interpCode[i];

    double highVar = fabs(high.getVal());
    double nominal = _nominal[i];
    double x = (high.getVal() < 0) ? -param.getVal() : +param.getVal();

    switch (icode)
    {

    case 0:
    {
      // piece-wise linear: only consider upper uncertainty as it should always be symmetric in XML workspace builder
      total *= nominal + x * highVar;
      break;
    }
    case 1:
    {
      // pice-wise log: only consider upper uncertainty as it should always be symmetric in XML workspace builder
      total *= logn(nominal, highVar, x);
      break;
    }
    case 4:
    {
      double boundary = _interpBoundary;
      double lowVar = fabs(low.getVal());

      if (x >= boundary)
      {
        total *= logn(nominal, highVar, x);
      }
      else if (x <= -boundary)
      {
        total *= logn(nominal, lowVar, x);
      }
      else if (x != 0)
      {
        total *= PolyInterpValue(i, x);
      }
      break;
    }
    default:
    {
      coutE(InputArguments) << "ResponseFunction::evaluate ERROR:  " << param.GetName()
                            << " with unknown interpolation code" << std::endl;
    }
    }
  }

  return total;
}

void ResponseFunction::printMultiline(std::ostream &os, Int_t contents,
                                      Bool_t verbose, TString indent) const
{
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  os << indent << "--- ResponseFunction ---" << std::endl;
  printResponseFunctions(os);
}

void ResponseFunction::printResponseFunctions(std::ostream &os) const
{
  for (int i = 0; i < _paramList.getSize(); i++)
  {
    const auto& param = static_cast<RooAbsReal&>(_paramList[i]);
    const auto& low = static_cast<RooAbsReal&>(_lowList[i]);
    const auto& high = static_cast<RooAbsReal&>(_highList[i]);
    os << std::setw(36) << param.GetName() << ": "
       << std::setw(7) << _nominal.at(i) << "  "
       << std::setw(7) << low.getVal() << "  "
       << std::setw(7) << high.getVal()
       << std::setw(7) << _interpCode.at(i)
       << std::endl;
  }
}

void ResponseFunction::cacheCoef(unsigned j) const
{
  // code for polynomial interpolation used when interpCode=4
  double boundary = _interpBoundary;
  double x0 = boundary;
    
  // location of the 6 coefficient for the j-th variable
  double *coeff = &_polCoeff[j * 6];
  double highVar = _nominal[j] + fabs((static_cast<RooAbsReal &>(_highList[j])).getVal());
  double lowVar = (_nominal[j] * _nominal[j]) / (_nominal[j] + fabs((static_cast<RooAbsReal &>(_lowList[j])).getVal()));

  // GHL: Swagato's suggestions
  double pow_up = std::pow(highVar / _nominal[j], x0);
  double pow_down = std::pow(lowVar / _nominal[j], x0);
  double logHi = std::log(highVar);
  double logLo = std::log(lowVar);
  double pow_up_log = highVar <= 0.0 ? 0.0 : pow_up * logHi;
  double pow_down_log = lowVar <= 0.0 ? 0.0 : -pow_down * logLo;
  double pow_up_log2 = highVar <= 0.0 ? 0.0 : pow_up_log * logHi;
  double pow_down_log2 = lowVar <= 0.0 ? 0.0 : -pow_down_log * logLo;

  double S0 = (pow_up + pow_down) / 2;
  double A0 = (pow_up - pow_down) / 2;
  double S1 = (pow_up_log + pow_down_log) / 2;
  double A1 = (pow_up_log - pow_down_log) / 2;
  double S2 = (pow_up_log2 + pow_down_log2) / 2;
  double A2 = (pow_up_log2 - pow_down_log2) / 2;

  // fcns+der+2nd_der are eq at bd

  // cache  coefficient of the polynomial
  coeff[0] = 1. / (8 * x0) * (15 * A0 - 7 * x0 * S1 + x0 * x0 * A2);
  coeff[1] = 1. / (8 * x0 * x0) * (-24 + 24 * S0 - 9 * x0 * A1 + x0 * x0 * S2);
  coeff[2] = 1. / (4 * pow(x0, 3)) * (-5 * A0 + 5 * x0 * S1 - x0 * x0 * A2);
  coeff[3] = 1. / (4 * pow(x0, 4)) * (12 - 12 * S0 + 7 * x0 * A1 - x0 * x0 * S2);
  coeff[4] = 1. / (8 * pow(x0, 5)) * (+3 * A0 - 3 * x0 * S1 + x0 * x0 * A2);
  coeff[5] = 1. / (8 * pow(x0, 6)) * (-8 + 8 * S0 - 5 * x0 * A1 + x0 * x0 * S2);
}

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,29,0)

#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFit/Detail/JSONInterface.h>
#include <RooFitHS3/JSONIO.h>
#include <RooWorkspace.h>

using namespace RooFit::Detail;

namespace {
  class ResponseFunctionStreamer : public RooFit::JSONIO::Exporter {
public:
  virtual const std::string& key() const override {
    static const std::string _key = "resonse";
    return _key;
  }
  
  virtual bool exportObject(RooJSONFactoryWSTool *, const RooAbsArg *func, JSONNode &elem) const override
  {
    const ResponseFunction *pip = static_cast<const ResponseFunction *>(func);
    elem["type"] << key();
    elem["interpolationCodes"].fill_seq(pip->interpolationCodes());
    elem["vars"].fill_seq(pip->parameters(), [](auto const &f) { return f->GetName(); });    
    elem["nom"].fill_seq(pip->nominal());
    elem["high"].fill_seq(pip->high(), [](auto const &f) { return f->GetName(); });
    elem["low"].fill_seq(pip->low(), [](auto const &f) { return f->GetName(); });    
    return true;
  }
};
} // namespace


namespace {
class ResponseFunctionFactory : public RooFit::JSONIO::Importer {
public:
   virtual bool importFunction(RooJSONFactoryWSTool *tool, const JSONNode &p) const override
   {
      std::string name(RooJSONFactoryWSTool::name(p));
      if (!p.has_child("vars")) {
         RooJSONFactoryWSTool::error("no vars of '" + name + "'");
      }
      if (!p.has_child("high")) {
         RooJSONFactoryWSTool::error("no high variations of '" + name + "'");
      }
      if (!p.has_child("low")) {
         RooJSONFactoryWSTool::error("no low variations of '" + name + "'");
      }
      if (!p.has_child("nom")) {
         RooJSONFactoryWSTool::error("no nominal variation of '" + name + "'");
      }

      std::vector<double> nom;
      nom << p["nom"];

      std::vector<int> codes;
      codes << p["interpolationCodes"];
      
      RooArgList vars;
      for (const auto &d : p["vars"].children()) {
         std::string objname(RooJSONFactoryWSTool::name(d));
         RooRealVar *obj = tool->request<RooRealVar>(objname, name);
         vars.add(*obj);
      }

      RooArgList high;
      for (const auto &d : p["high"].children()) {
         std::string objname(RooJSONFactoryWSTool::name(d));
         RooAbsReal *obj = tool->request<RooAbsReal>(objname, name);
         high.add(*obj);
      }

      RooArgList low;
      for (const auto &d : p["low"].children()) {
         std::string objname(RooJSONFactoryWSTool::name(d));
         RooAbsReal *obj = tool->request<RooAbsReal>(objname, name);
         low.add(*obj);
      }

      ResponseFunction pip(name.c_str(), name.c_str(), vars, nom, low, high, codes);
      
      tool->workspace()->import(pip, RooFit::RecycleConflictNodes(true), RooFit::Silence(true));
      return true;
   }
};
} // namespace


namespace {
  int register_serializations(){
    RooFit::JSONIO::registerImporter("response", std::make_unique<ResponseFunctionFactory>());
    RooFit::JSONIO::registerExporter(ResponseFunction::Class(), std::make_unique<ResponseFunctionStreamer>());
    return 1;
  }
  int _dummy = register_serializations();
}

#endif

