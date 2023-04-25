/** \class PiecewiseInterpolation
 * \ingroup HistFactory
 * The PiecewiseInterpolation is a class that can morph distributions into each other, which
 * is useful to estimate systematic uncertainties. Given a nominal distribution and one or
 * more altered or distorted ones, it computes a new shape depending on the value of the nuisance
 * parameters \f$ \alpha_i \f$:
 * \f[
 *   A = \sum_i \mathrm{Interpolate}(\mathrm{low}_i, \mathrm{nominal}, \mathrm{high}_i, \alpha_i).
 * \f]
 * If an \f$ \alpha_i \f$ is zero, the distribution is identical to the nominal distribution, at
 * \f$ \pm 1 \f$ it is identical to the up/down distribution for that specific \f$ i \f$.
 *
 * The class supports several interpolation methods, which can be selected for each parameter separately
 * using setInterpCode(). The default interpolation code is 4. This performs
 * - \f$ |\alpha | > 1 \f$: Linear extrapolation.
 * - \f$ |\alpha | < 1 \f$: Polynomial interpolation. A sixth-order polynomial is used. Its coefficients
 * are chosen such that function, first, and second derivative at \f$ \alpha \pm 1 \f$ match the values
 * that the extrapolation procedure uses.
 */

#include "RooStats/HistFactory/PiecewiseInterpolation.h"

#include "RooFit.h"

#include "Riostream.h"
#include "TBuffer.h"

#include "RooAbsReal.h"
#include "RooAbsSelfCachedReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooNumIntConfig.h"
#include "RooTrace.h"

#include <exception>
#include <math.h>

using namespace std;

ClassImp(PiecewiseInterpolation);
;


////////////////////////////////////////////////////////////////////////////////

PiecewiseInterpolation::PiecewiseInterpolation()
{
  _positiveDefinite=false;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Construct a new interpolation. The value of the function will be
/// \f[
///   A = \sum_i \mathrm{Interpolate}(\mathrm{low}_i, \mathrm{nominal}, \mathrm{high}_i).
/// \f]
/// \param name Name of the object.
/// \param title Title (for e.g. plotting)
/// \param nominal Nominal value of the function.
/// \param lowSet  Set of down variations.
/// \param highSet Set of up variations.
/// \param paramSet Parameters that control the interpolation.
/// \param takeOwnership If true, the PiecewiseInterpolation object will take ownership of the arguments in the low, high and parameter sets.
PiecewiseInterpolation::PiecewiseInterpolation(const char* name, const char* title, const RooAbsReal& nominal,
					       const RooArgList& lowSet, 
					       const RooArgList& highSet,
					       const RooArgList& paramSet,
					       Bool_t takeOwnership) :
  RooAbsSelfCachedReal(name, title, 0),
  _depList("!depList","list of observables", this),
  _nominal("!nominal","nominal value", this, (RooAbsReal&)nominal),
  _lowSet("!lowSet","low-side variation",this),
  _highSet("!highSet","high-side variation",this),
  _paramSet("!paramSet","high-side variation",this),
  _positiveDefinite(false)

{
  _depList.add(*(_nominal.arg().getObservables((((RooHistFunc&)_nominal.arg()).dataHist()))));
  // KC: check both sizes
  if (lowSet.getSize() != highSet.getSize()) {
    coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: input lists should be of equal length" << endl ;
    RooErrorHandler::softAbort() ;    
  }

  RooFIter inputIter1 = lowSet.fwdIterator() ;
  RooAbsArg* comp ;
  while((comp = inputIter1.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _lowSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
  }


  RooFIter inputIter2 = highSet.fwdIterator() ;
  while((comp = inputIter2.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _highSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
  }


  RooFIter inputIter3 = paramSet.fwdIterator() ;
  while((comp = inputIter3.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _paramSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
    _interpCode.push_back(0); // default code: linear interpolation
  }

  
  // Choose special integrator by default 
  specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;
  TRACE_CREATE
  // this->Print("V");
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

PiecewiseInterpolation::PiecewiseInterpolation(const PiecewiseInterpolation& other, const char* name) :
  RooAbsSelfCachedReal(other, name),
  _depList("!depList",this,other._depList),
  _nominal("!nominal",this,other._nominal),
  _lowSet("!lowSet",this,other._lowSet),
  _highSet("!highSet",this,other._highSet),
  _paramSet("!paramSet",this,other._paramSet),
  _positiveDefinite(other._positiveDefinite),
  _interpCode(other._interpCode)
{
  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

PiecewiseInterpolation::~PiecewiseInterpolation() 
{
  TRACE_DESTROY
}



RooArgSet* PiecewiseInterpolation::actualObservables(const RooArgSet &nset) const
{
  RooArgSet *myDeps = new RooArgSet;
  for (auto comp : _depList) {
    auto dep = static_cast<RooAbsArg*>(comp);
    myDeps->add(*dep);
  }
  return myDeps;
}


RooArgSet* PiecewiseInterpolation::actualParameters(const RooArgSet &nset) const
{
  RooArgSet *myPars=new RooArgSet;
  for (auto comp : _paramSet) {
    auto par = static_cast<RooAbsArg*>(comp);
    myPars->add(*par);
  }
  myPars->remove(nset, kTRUE, kTRUE);
  return myPars;
}


Double_t PiecewiseInterpolation::getValV(const RooArgSet* nset) const
{
  return RooAbsCachedReal::getValV(((RooArgSet*)&_depList));
}


void PiecewiseInterpolation::fillCacheObject(RooAbsCachedReal::FuncCacheElem& cache) const
{
  RooDataHist& cacheHist = *(cache.hist()) ;
  RooDataHist& nomHist=((RooHistFunc&)(_nominal.arg())).dataHist();

  // Iterator over all bins of RooDataHist and fill weights
  for (Int_t k=0 ; k<cacheHist.numEntries() ; k++) {
    cacheHist.get(k) ;
    nomHist.get(k);

    Double_t nominal = nomHist.weight();
    Double_t sum(nominal) ;

  //  return sum;
    RooAbsReal* param ;
    RooHistFunc* high ;
    RooHistFunc* low ;
    int i=-1;

    RooFIter lowIter(_lowSet.fwdIterator()) ;
    RooFIter highIter(_highSet.fwdIterator()) ;
    RooFIter paramIter(_paramSet.fwdIterator()) ;

    while((param=(RooAbsReal*)paramIter.next())) {
      ++i;
      low = (RooHistFunc*)(lowIter.next()) ;
      high = (RooHistFunc*)(highIter.next()) ;
      low->dataHist().get(k);
      high->dataHist().get(k);
      Double_t highVal=high->dataHist().weight();
      Double_t lowVal=low->dataHist().weight();
      if ((lowVal==nominal)*(highVal==nominal))
      {
        continue;
      }

      Int_t icode = _interpCode[i] ;

      switch(icode) {
      case 0: {
        // piece-wise linear
        if(param->getVal()>0)
    sum +=  param->getVal()*(highVal - nominal );
        else
    sum += param->getVal()*(nominal - lowVal);
        break ;
      }
      case 1: {
        // pice-wise log
        if(param->getVal()>=0)
    sum *= pow(highVal/nominal, +param->getVal());
        else
    sum *= pow(lowVal/nominal,  -param->getVal());
        break ;
      }
      case 2: {
        // parabolic with linear
        double a = 0.5*(highVal+lowVal)-nominal;
        double b = 0.5*(highVal-lowVal);
        double c = 0;
        double x = param->getVal();
        if(x*x <= 1)
        {
          sum += a*x*x+b*x;
        }
        else
        {
          int sgn = 1-2*std::signbit(x);
          sum += (b+2*sgn*a)*x-a;
        }
        /* can I do this with no branch?
        if(param->getVal()>1 ){
    sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
        } else if(param->getVal()<-1 ) {
    sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
        } else {
    sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
        }
        */
        break ;
      }
      case 3: {
        //parabolic version of log-normal
        double a = 0.5*(highVal+lowVal)-nominal;
        double b = 0.5*(highVal-lowVal);
        double c = 0;
        if(param->getVal()>1 ){
    sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
        } else if(param->getVal()<-1 ) {
    sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
        } else {
    sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
        }
        break ;
      }
      case 4: {

        // WVE ****************************************************************
        // WVE *** THIS CODE IS CRITICAL TO HISTFACTORY FIT CPU PERFORMANCE ***
        // WVE *** Do not modify unless you know what you are doing...      ***
        // WVE ****************************************************************

        double x  = param->getVal();
        if (x>1) {
    sum += x*(highVal - nominal );
        } else if (x<-1) {
    sum += x*(nominal - lowVal);
        } else {
    double eps_plus = highVal - nominal;
    double eps_minus = nominal - lowVal;
    double S = 0.5 * (eps_plus + eps_minus);
    double A = 0.0625 * (eps_plus - eps_minus);
    double val = nominal + x * (S + x * A * ( 15 + x * x * (-10 + x * x * 3  ) ) );

    if (val < 0) val = 0;
    sum += val-nominal;
        }
        break ;

        // WVE ****************************************************************
      }
      case 5: {

        double x0 = 1.0;//boundary;
        double x  = param->getVal();

        if (x > x0 || x < -x0)
        {
    if(x>0)
      sum += x*(highVal - nominal );
    else
      sum += x*(nominal - lowVal);
        }
        else if (nominal != 0)
        {
    double eps_plus = highVal - nominal;
    double eps_minus = nominal - lowVal;
    double S = (eps_plus + eps_minus)/2;
    double A = (eps_plus - eps_minus)/2;

    //fcns+der are eq at bd
    double a = S;
    double b = 3*A/(2*x0);
    //double c = 0;
    double d = -A/(2*x0*x0*x0);

    double val = nominal + a*x + b*pow(x, 2) + 0/*c*pow(x, 3)*/ + d*pow(x, 4);
    if (val < 0) val = 0;

    //cout << "Using interp code 5, val = " << val << endl;

    sum += val-nominal;
        }
        break ;
      }
      default: {
        coutE(InputArguments) << "PiecewiseInterpolation::evaluate ERROR:  " << param->GetName() 
            << " with unknown interpolation code" << icode << endl ;
        break ;
      }
      }
    }

    if(_positiveDefinite && (sum<0)){
      sum = 1e-6;
      sum = 0;
      //     cout <<"sum < 0 forcing  positive definite"<<endl;
      //     int code = 1;
      //     RooArgSet* myset = new RooArgSet();
      //     cout << "integral = " << analyticalIntegralWN(code, myset) << endl;
    } else if(sum<0){
       cxcoutD(Tracing) <<"PiecewiseInterpolation::evaluate -  sum < 0, not forcing positive definite"<<endl;
    }
    cacheHist.set(sum);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate and return current value of self

Double_t PiecewiseInterpolation::evaluate() const 
{
  ///////////////////
  Double_t nominal = _nominal;
  Double_t sum(nominal) ;

  for (unsigned int i=0; i < _paramSet.size(); ++i) {
    auto param = static_cast<RooAbsReal*>(_paramSet.at(i));
    auto low   = static_cast<RooAbsReal*>(_lowSet.at(i));
    auto high  = static_cast<RooAbsReal*>(_highSet.at(i));
    Int_t icode = _interpCode[i] ;

    switch(icode) {
    case 0: {
      // piece-wise linear
      if(param->getVal()>0)
	sum +=  param->getVal()*(high->getVal() - nominal );
      else
	sum += param->getVal()*(nominal - low->getVal());
      break ;
    }
    case 1: {
      // pice-wise log
      if(param->getVal()>=0)
	sum *= pow(high->getVal()/nominal, +param->getVal());
      else
	sum *= pow(low->getVal()/nominal,  -param->getVal());
      break ;
    }
    case 2: {
      // parabolic with linear
      double a = 0.5*(high->getVal()+low->getVal())-nominal;
      double b = 0.5*(high->getVal()-low->getVal());
      double c = 0;
      if(param->getVal()>1 ){
	sum += (2*a+b)*(param->getVal()-1)+high->getVal()-nominal;
      } else if(param->getVal()<-1 ) {
	sum += -1*(2*a-b)*(param->getVal()+1)+low->getVal()-nominal;
      } else {
	sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }
    case 3: {
      //parabolic version of log-normal
      double a = 0.5*(high->getVal()+low->getVal())-nominal;
      double b = 0.5*(high->getVal()-low->getVal());
      double c = 0;
      if(param->getVal()>1 ){
	sum += (2*a+b)*(param->getVal()-1)+high->getVal()-nominal;
      } else if(param->getVal()<-1 ) {
	sum += -1*(2*a-b)*(param->getVal()+1)+low->getVal()-nominal;
      } else {
	sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }
    case 4: {
      
      // WVE ****************************************************************
      // WVE *** THIS CODE IS CRITICAL TO HISTFACTORY FIT CPU PERFORMANCE ***
      // WVE *** Do not modify unless you know what you are doing...      ***
      // WVE ****************************************************************

      double x  = param->getVal();      
      if (x>1) {
	sum += x*(high->getVal() - nominal );
      } else if (x<-1) {
	sum += x*(nominal - low->getVal());
      } else {
	double eps_plus = high->getVal() - nominal;
	double eps_minus = nominal - low->getVal();
	double S = 0.5 * (eps_plus + eps_minus);
	double A = 0.0625 * (eps_plus - eps_minus);
	
	//fcns+der+2nd_der are eq at bd

        double val = nominal + x * (S + x * A * ( 15 + x * x * (-10 + x * x * 3  ) ) ); 


	if (val < 0) val = 0;
	sum += val-nominal;
      }
      break ;

      // WVE ****************************************************************
    }
    case 5: {
      
      double x0 = 1.0;//boundary;
      double x  = param->getVal();

      if (x > x0 || x < -x0)
      {
	if(x>0)
	  sum += x*(high->getVal() - nominal );
	else
	  sum += x*(nominal - low->getVal());
      }
      else if (nominal != 0)
      {
	double eps_plus = high->getVal() - nominal;
	double eps_minus = nominal - low->getVal();
	double S = (eps_plus + eps_minus)/2;
	double A = (eps_plus - eps_minus)/2;

	//fcns+der are eq at bd
	double a = S;
	double b = 3*A/(2*x0);
	//double c = 0;
	double d = -A/(2*x0*x0*x0);

	double val = nominal + a*x + b*pow(x, 2) + 0/*c*pow(x, 3)*/ + d*pow(x, 4);
	if (val < 0) val = 0;

	//cout << "Using interp code 5, val = " << val << endl;

	sum += val-nominal;
      }
      break ;
    }
    default: {
      coutE(InputArguments) << "PiecewiseInterpolation::evaluate ERROR:  " << param->GetName() 
			    << " with unknown interpolation code" << icode << endl ;
      break ;
    }
    }
  }
  
  if(_positiveDefinite && (sum<0)){
    sum = 1e-6;
    sum = 0;
    //     cout <<"sum < 0 forcing  positive definite"<<endl;
     //     int code = 1;
     //     RooArgSet* myset = new RooArgSet();
     //     cout << "integral = " << analyticalIntegralWN(code, myset) << endl;
  } else if(sum<0){
     cxcoutD(Tracing) <<"PiecewiseInterpolation::evaluate -  sum < 0, not forcing positive definite"<<endl;
  }
  return sum;

}

////////////////////////////////////////////////////////////////////////////////

Bool_t PiecewiseInterpolation::setBinIntegrator(RooArgSet& allVars) 
{
  if(allVars.getSize()==1){
    RooAbsReal* temp = const_cast<PiecewiseInterpolation*>(this);
    temp->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator")  ;
    int nbins = ((RooRealVar*) allVars.first())->numBins();
    temp->specialIntegratorConfig(kTRUE)->getConfigSection("RooBinIntegrator").setRealValue("numBins",nbins);
    return true;
  }else{
    cout << "Currently BinIntegrator only knows how to deal with 1-d "<<endl;
    return false;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise that all integrals can be handled internally.


////////////////////////////////////////////////////////////////////////////////
/// Implement analytical integrations by doing appropriate weighting from  component integrals
/// functions to integrators of components

#if 0
Double_t PiecewiseInterpolation::analyticalIntegralWN(Int_t code, const RooArgSet* /*normSet2*/,const char* /*rangeName*/) const 
{
  /*
  cout <<"Enter analytic Integral"<<endl;
  printDirty(true);
  //  _nominal.arg().setDirtyInhibit(kTRUE) ;
  _nominal.arg().setShapeDirty() ;
  RooAbsReal* temp ;
  RooFIter lowIter(_lowSet.fwdIterator()) ;
  while((temp=(RooAbsReal*)lowIter.next())) {
    //    temp->setDirtyInhibit(kTRUE) ;
    temp->setShapeDirty() ;
  }
  RooFIter highIter(_highSet.fwdIterator()) ;
  while((temp=(RooAbsReal*)highIter.next())) {
    //    temp->setDirtyInhibit(kTRUE) ;
    temp->setShapeDirty() ;
  }
  */

  /*
  RooAbsArg::setDirtyInhibit(kTRUE);
  printDirty(true);
  cout <<"done setting dirty inhibit = true"<<endl;

  // old integral, only works for linear and not positive definite
  CacheElem* cache = (CacheElem*) _normIntMgr.getObjByIndex(code-1) ;

  
 std::unique_ptr<RooArgSet> vars2( getParameters(RooArgSet()) );
 std::unique_ptr<RooArgSet> iset(  _normIntMgr.nameSet2ByIndex(code-1)->select(*vars2) );            
 cout <<"iset = "<<endl;
 iset->Print("v");

  double sum = 0;
  RooArgSet* vars = getVariables();
  vars->remove(_paramSet);
  _paramSet.Print("v");
  vars->Print("v");
  if(vars->getSize()==1){
    RooRealVar* obs = (RooRealVar*) vars->first();
    for(int i=0; i<obs->numBins(); ++i){
      obs->setVal( obs->getMin() + (.5+i)*(obs->getMax()-obs->getMin())/obs->numBins());
      sum+=evaluate()*(obs->getMax()-obs->getMin())/obs->numBins();
      cout << "obs = " << obs->getVal() << " sum = " << sum << endl;
    }
  } else{
    cout <<"only know how to deal with 1 observable right now"<<endl;
  }
  */

  /*
  _nominal.arg().setDirtyInhibit(kFALSE) ;
  RooFIter lowIter2(_lowSet.fwdIterator()) ;
  while((temp=(RooAbsReal*)lowIter2.next())) {
    temp->setDirtyInhibit(kFALSE) ;
  }
  RooFIter highIter2(_highSet.fwdIterator()) ;
  while((temp=(RooAbsReal*)highIter2.next())) {
    temp->setDirtyInhibit(kFALSE) ;
  }
  */
  
  /*
  RooAbsArg::setDirtyInhibit(kFALSE);
  printDirty(true);
  cout <<"done"<<endl;
  cout << "sum = " <<sum<<endl;
  //return sum;
  */  

  // old integral, only works for linear and not positive definite
  CacheElem* cache = (CacheElem*) _normIntMgr.getObjByIndex(code-1) ;
  if( cache==NULL ) {
    std::cout << "Error: Cache Element is NULL" << std::endl;
    throw std::exception();
  }

  // old integral, only works for linear and not positive definite
  RooFIter funcIntIter = cache->_funcIntList.fwdIterator() ;
  RooFIter lowIntIter = cache->_lowIntList.fwdIterator() ;
  RooFIter highIntIter = cache->_highIntList.fwdIterator() ;
  RooAbsReal *funcInt(0), *low(0), *high(0), *param(0) ;
  Double_t value(0) ;
  Double_t nominal(0);

  // get nominal 
  int i=0;
  while(( funcInt = (RooAbsReal*)funcIntIter.next())) {
    value += funcInt->getVal() ;
    nominal = value;
    i++;
  }
  if(i==0 || i>1) { cout << "problem, wrong number of nominal functions"<<endl; }

  // now get low/high variations
  i = 0;
  RooFIter paramIter(_paramSet.fwdIterator()) ;

  // KC: old interp code with new iterator
  while( (param=(RooAbsReal*)paramIter.next()) ) {
    low = (RooAbsReal*)lowIntIter.next() ;
    high = (RooAbsReal*)highIntIter.next() ;
    
    if(param->getVal()>0) {
      value += param->getVal()*(high->getVal() - nominal );
    } else {
      value += param->getVal()*(nominal - low->getVal());
    }
    ++i;
  }

  /* // MB : old bit of interpolation code
  while( (param=(RooAbsReal*)_paramIter->Next()) ) {
    low = (RooAbsReal*)lowIntIter->Next() ;
    high = (RooAbsReal*)highIntIter->Next() ;
    
    if(param->getVal()>0) {
      value += param->getVal()*(high->getVal() - nominal );
    } else {
      value += param->getVal()*(nominal - low->getVal());
    }
    ++i;
  }
  */

  /* KC: the code below is wrong.  Can't pull out a constant change to a non-linear shape deformation.
  while( (param=(RooAbsReal*)paramIter.next()) ) {
    low = (RooAbsReal*)lowIntIter.next() ;
    high = (RooAbsReal*)highIntIter.next() ;

    if(_interpCode.empty() || _interpCode.at(i)==0){
      // piece-wise linear
      if(param->getVal()>0)
	value +=  param->getVal()*(high->getVal() - nominal );
      else
	value += param->getVal()*(nominal - low->getVal());
    } else if(_interpCode.at(i)==1){
      // pice-wise log
      if(param->getVal()>=0)
	value *= pow(high->getVal()/nominal, +param->getVal());
      else
	value *= pow(low->getVal()/nominal,  -param->getVal());
    } else if(_interpCode.at(i)==2){
      // parabolic with linear
      double a = 0.5*(high->getVal()+low->getVal())-nominal;
      double b = 0.5*(high->getVal()-low->getVal());
      double c = 0;
      if(param->getVal()>1 ){
	value += (2*a+b)*(param->getVal()-1)+high->getVal()-nominal;
      } else if(param->getVal()<-1 ) {
	value += -1*(2*a-b)*(param->getVal()+1)+low->getVal()-nominal;
      } else {
	value +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
    } else if(_interpCode.at(i)==3){
      //parabolic version of log-normal
      double a = 0.5*(high->getVal()+low->getVal())-nominal;
      double b = 0.5*(high->getVal()-low->getVal());
      double c = 0;
      if(param->getVal()>1 ){
	value += (2*a+b)*(param->getVal()-1)+high->getVal()-nominal;
      } else if(param->getVal()<-1 ) {
	value += -1*(2*a-b)*(param->getVal()+1)+low->getVal()-nominal;
      } else {
	value +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
	
    } else {
      coutE(InputArguments) << "PiecewiseInterpolation::analyticalIntegralWN ERROR:  " << param->GetName() 
			    << " with unknown interpolation code" << endl ;
    }
    ++i;
  }
  */

  //  cout << "value = " << value <<endl;
  return value;
}
#endif


////////////////////////////////////////////////////////////////////////////////

void PiecewiseInterpolation::setInterpCode(RooAbsReal& param, int code){
  int index = _paramSet.index(&param);
  if(index<0){
      coutE(InputArguments) << "PiecewiseInterpolation::setInterpCode ERROR:  " << param.GetName() 
			    << " is not in list" << endl ;
  } else {
      coutW(InputArguments) << "PiecewiseInterpolation::setInterpCode :  " << param.GetName() 
			    << " is now " << code << endl ;
    _interpCode.at(index) = code;
  }
}


////////////////////////////////////////////////////////////////////////////////

void PiecewiseInterpolation::setAllInterpCodes(int code){
  for(unsigned int i=0; i<_interpCode.size(); ++i){
    _interpCode.at(i) = code;
  }
}


////////////////////////////////////////////////////////////////////////////////

void PiecewiseInterpolation::printAllInterpCodes(){
  for(unsigned int i=0; i<_interpCode.size(); ++i){
    coutI(InputArguments) <<"interp code for " << _paramSet.at(i)->GetName() << " = " << _interpCode.at(i) <<endl;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// WVE note: assumes nominal and alternates have identical structure, must add explicit check

std::list<Double_t>* PiecewiseInterpolation::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  return _nominal.arg().binBoundaries(obs,xlo,xhi) ;  
}


////////////////////////////////////////////////////////////////////////////////
/// WVE note: assumes nominal and alternates have identical structure, must add explicit check

Bool_t PiecewiseInterpolation::isBinnedDistribution(const RooArgSet& obs) const 
{
  return _nominal.arg().isBinnedDistribution(obs) ;
}



////////////////////////////////////////////////////////////////////////////////

std::list<Double_t>* PiecewiseInterpolation::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  return _nominal.arg().plotSamplingHint(obs,xlo,xhi) ;  
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class PiecewiseInterpolation.

void PiecewiseInterpolation::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PiecewiseInterpolation::Class(),this);
      specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;      
      if (_interpCode.empty()) _interpCode.resize(_paramSet.getSize());
   } else {
      R__b.WriteClassBuffer(PiecewiseInterpolation::Class(),this);
   }
}


/*
////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a PiecewiseInterpolation to more intuitively reflect the contents of the
/// product operator construction

void PiecewiseInterpolation::printMetaArgs(ostream& os) const 
{
  _lowIter->Reset() ;
  if (_highIter) {
    _highIter->Reset() ;
  }

  Bool_t first(kTRUE) ;
    
  RooAbsArg* arg1, *arg2 ;
  if (_highSet.getSize()!=0) { 

    while((arg1=(RooAbsArg*)_lowIter->Next())) {
      if (!first) {
	os << " + " ;
      } else {
	first = kFALSE ;
      }
      arg2=(RooAbsArg*)_highIter->Next() ;
      os << arg1->GetName() << " * " << arg2->GetName() ;
    }

  } else {
    
    while((arg1=(RooAbsArg*)_lowIter->Next())) {
      if (!first) {
	os << " + " ;
      } else {
	first = kFALSE ;
      }
      os << arg1->GetName() ; 
    }  

  }

  os << " " ;    
}

*/
