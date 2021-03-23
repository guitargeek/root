#include "RooFitDriver.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"

#include "RooBatchCompute.h"

#include <numeric>
#include <stdexcept>
#include <vector>

ClassImp(RooFitDriver);

RooFitDriver::RooFitDriver(const char *name, const char *title, RooAbsReal& absReal, RooAbsData& data, RooArgSet const& normSet)
  : RooAbsReal(name,title), _absReal{&absReal}, _data{&data}, _normSet{normSet}
{}

RooFitDriver::RooFitDriver(const RooFitDriver& other, const char* name)
  : RooAbsReal(other, name), _absReal{other._absReal}, _data{other._data}, _normSet{other._normSet}
{}

double RooFitDriver::evaluate() const {
  throw std::runtime_error("RooFitDriver::evaluate was called directly which should not happen!");
}

double RooFitDriver::getValV(const RooArgSet* normSet) const {

  if(normSet) {
    throw std::runtime_error("RooFitDriver::getValV should always be passed a nullptr normSet!");
  }

  // smart run context clearing logic should go here
  _runContext.clear();

  // get the data into the run context
  _data->getBatches(_runContext, 0, _data->numEntries());

  _absReal->getValues(_runContext, &_normSet);

  // this is the only time where data should be copied from the RunContext
  auto values = _runContext[_absReal];

  if(values.size() != 1) {
    throw std::runtime_error("RooFitDriver was used with a RooAbsReal that is not a scalar!");
  }

  return values[0];
}

RooSpan<double> RooFitDriver::evaluateSpan(RooBatchCompute::RunContext&, const RooArgSet*) const {
  throw std::runtime_error("RooFitDriver::evaluateSpan was called which should not happen!");
}

RooArgSet* RooFitDriver::getParameters(const RooArgSet*, Bool_t stripDisconnected) const {
  return _absReal->getParameters(*_data, stripDisconnected);
}
