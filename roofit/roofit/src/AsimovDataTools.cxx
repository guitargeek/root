/*
 * Project: RooFit
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include <RooFit/AsimovDataTools.h>

#include <RooCategory.h>
#include <RooGaussian.h>
#include <RooPoisson.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>

#include <TMath.h>

namespace {

////////////////////////////////////////////////////////////////////////////////
/// set observed value to the expected one
/// works for Gaussian, Poisson or LogNormal
/// assumes mean parameter value is the argument not constant and not depending on observables
/// (if more than two arguments are not constant will use first one but print a warning !)
/// need to iterate on the components of the Poisson to get n and nu (nu can be a RooAbsReal)
/// (code from G. Petrucciani and extended by L.M.)

bool setObsToExpected(RooAbsPdf &pdf, const RooArgSet &obs, int printLevel)
{
   RooRealVar *myobs = 0;
   RooAbsReal *myexp = 0;
   const char *pdfName = pdf.IsA()->GetName();
   RooFIter iter(pdf.serverMIterator());
   for (RooAbsArg *a = iter.next(); a != 0; a = iter.next()) {
      if (obs.contains(*a)) {
         if (myobs != 0) {
            oocoutF((TObject *)0, Generation)
               << "AsymptoticCalculator::SetObsExpected( " << pdfName << " ) : Has two observables ?? " << std::endl;
            return false;
         }
         myobs = dynamic_cast<RooRealVar *>(a);
         if (myobs == 0) {
            oocoutF((TObject *)0, Generation) << "AsymptoticCalculator::SetObsExpected( " << pdfName
                                              << " ) : Observable is not a RooRealVar??" << std::endl;
            return false;
         }
      } else {
         if (!a->isConstant()) {
            if (myexp != 0) {
               oocoutE((TObject *)0, Generation) << "AsymptoticCalculator::SetObsExpected( " << pdfName
                                                 << " ) : Has two non-const arguments  " << std::endl;
               return false;
            }
            myexp = dynamic_cast<RooAbsReal *>(a);
            if (myexp == 0) {
               oocoutF((TObject *)0, Generation) << "AsymptoticCalculator::SetObsExpected( " << pdfName
                                                 << " ) : Expected is not a RooAbsReal??" << std::endl;
               return false;
            }
         }
      }
   }
   if (myobs == 0) {
      oocoutF((TObject *)0, Generation) << "AsymptoticCalculator::SetObsExpected( " << pdfName << " ) : No observable?"
                                        << std::endl;
      return false;
   }
   if (myexp == 0) {
      oocoutF((TObject *)0, Generation) << "AsymptoticCalculator::SetObsExpected( " << pdfName << " ) : No observable?"
                                        << std::endl;
      return false;
   }

   myobs->setVal(myexp->getVal());

   if (printLevel > 2) {
      std::cout << "setObsToExpected : setting " << myobs->GetName() << " to expected value " << myexp->getVal()
                << " of " << myexp->GetName() << std::endl;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Inpspect a product pdf to find all the Poisson or Gaussian parts to set the observed
/// values to expected ones.

bool setObsToExpected(RooProdPdf &prod, const RooArgSet &obs, int printLevel)
{
   RooLinkedListIter iter(prod.pdfList().iterator());
   bool ret = true;
   for (RooAbsArg *a = (RooAbsArg *)iter.Next(); a != 0; a = (RooAbsArg *)iter.Next()) {
      if (!a->dependsOn(obs))
         continue;
      RooPoisson *pois = 0;
      RooGaussian *gaus = 0;
      if ((pois = dynamic_cast<RooPoisson *>(a)) != 0) {
         ret &= setObsToExpected(*pois, obs, printLevel);
         pois->setNoRounding(true); // needed since expected value is not an integer
      } else if ((gaus = dynamic_cast<RooGaussian *>(a)) != 0) {
         ret &= setObsToExpected(*gaus, obs, printLevel);
      } else {
         // should try to add also lognormal case ?
         RooProdPdf *subprod = dynamic_cast<RooProdPdf *>(a);
         if (subprod)
            ret &= setObsToExpected(*subprod, obs, printLevel);
         else {
            oocoutE((TObject *)0, InputArguments)
               << "Illegal term in counting model: "
               << "the PDF " << a->GetName() << " depends on the observables, but is not a Poisson, Gaussian or Product"
               << std::endl;
            return false;
         }
      }
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// fill bins by looping recursively on observables

void fillBins(const RooAbsPdf &pdf, const RooArgList &obs, RooAbsData &data, int &index, double &binVolume, int &ibin,
              int printLevel)
{

   bool debug = (printLevel >= 2);

   RooRealVar *v = dynamic_cast<RooRealVar *>(&(obs[index]));
   if (!v)
      return;

   RooArgSet obstmp(obs);
   double expectedEvents = pdf.expectedEvents(obstmp);
   // if (debug)  {
   //    std::cout << "expected events = " << expectedEvents << std::endl;
   // }

   if (debug)
      std::cout << "looping on observable " << v->GetName() << std::endl;
   for (int i = 0; i < v->getBins(); ++i) {
      v->setBin(i);
      if (index < obs.getSize() - 1) {
         index++; // increase index
         double prevBinVolume = binVolume;
         binVolume *= v->getBinWidth(i); // increase bin volume
         fillBins(pdf, obs, data, index, binVolume, ibin, printLevel);
         index--;                   // decrease index
         binVolume = prevBinVolume; // decrease also bin volume
      } else {

         // this is now a new bin - compute the pdf in this bin
         double totBinVolume = binVolume * v->getBinWidth(i);
         double fval = pdf.getVal(&obstmp) * totBinVolume;

         // if (debug) std::cout << "pdf value in the bin " << fval << " bin volume = " << totBinVolume << "   " <<
         // fval*expectedEvents << std::endl;
         if (fval * expectedEvents <= 0) {
            if (fval * expectedEvents < 0) {
               oocoutW(static_cast<TObject *>(nullptr), InputArguments)
                  << "AsymptoticCalculator::" << __func__
                  << "(): Detected a bin with negative expected events! Please check your inputs." << std::endl;
            } else {
               oocoutW(static_cast<TObject *>(nullptr), InputArguments)
                  << "AsymptoticCalculator::" << __func__ << "(): Detected a bin with zero expected events- skip it"
                  << std::endl;
            }
         }
         // have a cut off for overflows ??
         else
            data.add(obs, fval * expectedEvents);

         if (debug) {
            std::cout << "bin " << ibin << "\t";
            for (int j = 0; j < obs.getSize(); ++j) {
               std::cout << "  " << ((RooRealVar &)obs[j]).getVal();
            }
            std::cout << " w = " << fval * expectedEvents;
            std::cout << std::endl;
         }
         // RooArgSet xxx(obs);
         // h3->Fill(((RooRealVar&) obs[0]).getVal(), ((RooRealVar&) obs[1]).getVal(), ((RooRealVar&) obs[2]).getVal() ,
         //          pdf->getVal(&xxx) );
         ibin++;
      }
   }
   // reset bin values
   if (debug)
      std::cout << "ending loop on .. " << v->GetName() << std::endl;

   v->setBin(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate counting Asimov data for the case when the pdf cannot be extended.
/// This function assumes that the pdf is a RooPoisson or can be decomposed in a product of RooPoisson,
/// or is a RooGaussian. Otherwise, we cannot know how to make the Asimov data sets.

RooDataSet *generateCountingAsimov(RooAbsPdf &pdf, const RooArgSet &observables, const RooRealVar &,
                                   RooCategory *channelCat, int printLevel)
{
   RooArgSet obs(observables);
   RooProdPdf *prod = dynamic_cast<RooProdPdf *>(&pdf);
   RooPoisson *pois = 0;
   RooGaussian *gaus = 0;

   if (printLevel > 1)
      std::cout << "generate counting Asimov data for pdf of type " << pdf.IsA()->GetName() << std::endl;

   bool r = false;
   if (prod != 0) {
      r = setObsToExpected(*prod, observables, printLevel);
   } else if ((pois = dynamic_cast<RooPoisson *>(&pdf)) != 0) {
      r = setObsToExpected(*pois, observables, printLevel);
      // we need in this case to set Poisson to real values
      pois->setNoRounding(true);
   } else if ((gaus = dynamic_cast<RooGaussian *>(&pdf)) != 0) {
      r = setObsToExpected(*gaus, observables, printLevel);
   } else {
      oocoutE((TObject *)0, InputArguments)
         << "A counting model pdf must be either a RooProdPdf or a RooPoisson or a RooGaussian" << std::endl;
   }
   if (!r)
      return 0;
   int icat = 0;
   if (channelCat) {
      icat = channelCat->getCurrentIndex();
   }

   RooDataSet *ret = new RooDataSet(std::string("CountingAsimovData") + std::to_string(icat),
                                    std::string("CountingAsimovData") + std::to_string(icat), obs);
   ret->add(obs);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the asimov data set for an observable of a pdf.
/// It generates binned data following the binning of the observables.
// TODO: (possibility to change number of bins)
// TODO: implement integration over bin content

RooDataSet *generateAsimovSinglePdf(const RooAbsPdf &pdf, const RooArgSet &allobs, const RooRealVar &weightVar,
                                    RooCategory *channelCat, int printLevel)
{

   // Get observables defined by the pdf associated with this state
   std::unique_ptr<RooArgSet> obs(pdf.getObservables(allobs));

   // if pdf cannot be extended assume is then a counting experiment
   if (!pdf.canBeExtended())
      return generateCountingAsimov(const_cast<RooAbsPdf &>(pdf), *obs, weightVar, channelCat, printLevel);

   RooArgSet obsAndWeight(*obs);
   obsAndWeight.add(weightVar);

   RooDataSet *asimovData = 0;
   if (channelCat) {
      int icat = channelCat->getCurrentIndex();
      asimovData = new RooDataSet(std::string("AsimovData") + std::to_string(icat),
                                  std::string("combAsimovData") + std::to_string(icat),
                                  RooArgSet(obsAndWeight, *channelCat), RooFit::WeightVar(weightVar));
   } else
      asimovData = new RooDataSet("AsimovData", "AsimovData", RooArgSet(obsAndWeight), RooFit::WeightVar(weightVar));

   // This works only for 1D observables
   // RooRealVar* thisObs = ((RooRealVar*)obstmp->first());

   RooArgList obsList(*obs);

   // loop on observables and on the bins
   if (printLevel >= 2) {
      std::cout << "Generating Asimov data for pdf " << pdf.GetName() << std::endl;
      std::cout << "list of observables  " << std::endl;
      obsList.Print();
   }

   int obsIndex = 0;
   double binVolume = 1;
   int nbins = 0;
   fillBins(pdf, obsList, *asimovData, obsIndex, binVolume, nbins, printLevel);
   if (printLevel >= 2)
      std::cout << "filled from " << pdf.GetName() << "   " << nbins << " nbins "
                << " volume is " << binVolume << std::endl;

   // for (int iobs = 0; iobs < obsList.getSize(); ++iobs) {
   //    RooRealVar * thisObs = dynamic_cast<RooRealVar*> &obsList[i];
   //    if (thisObs == 0) continue;
   //    // loop on the bin contents
   //    for(int  ibin=0; ibin<thisObs->numBins(); ++ibin){
   //       thisObs->setBin(ibin);

   //   thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
   //   if (thisNorm*expectedEvents <= 0)
   //   {
   //     std::cout << "WARNING::Detected bin with zero expected events! Please check your inputs." << std::endl;
   //   }
   //   // have a cut off for overflows ??
   //   obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
   // }

   if (printLevel >= 1) {
      asimovData->Print();
      // cout <<"sum entries "<< asimovData->sumEntries()<<endl;
   }
   if (TMath::IsNaN(asimovData->sumEntries())) {
      std::cout << "sum entries is nan" << std::endl;
      assert(0);
      delete asimovData;
      asimovData = 0;
   }

   return asimovData;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Generate the Asimov data for the observables (not the global ones)
/// need to deal with the case of a sim pdf.

std::unique_ptr<RooDataSet>
RooFit::generateAsimovDataset(RooAbsPdf const &pdf, RooArgSet const &observables, int printLevel)
{

   auto weightVar = std::make_unique<RooRealVar>("binWeightAsimov", "binWeightAsimov", 1, 0, 1.E30);

   if (printLevel > 1) {
      std::cout << " Generate Asimov data for observables" << std::endl;
   }
   auto *simPdf = dynamic_cast<const RooSimultaneous *>(&pdf);
   if (!simPdf) {
      // generate data for non sim pdf
      return std::unique_ptr<RooDataSet>{generateAsimovSinglePdf(pdf, observables, *weightVar, 0, printLevel)};
   }

   std::map<std::string, RooDataSet *> asimovDataMap;

   // look at category of simpdf
   RooCategory &channelCat = const_cast<RooCategory &>(dynamic_cast<const RooCategory &>(simPdf->indexCat()));
   int nrIndices = channelCat.numTypes();
   if (nrIndices == 0) {
      oocoutW((TObject *)0, Generation) << "Simultaneous pdf does not contain any categories." << std::endl;
   }
   for (int i = 0; i < nrIndices; i++) {
      channelCat.setIndex(i);
      // Get pdf associated with state from simpdf
      RooAbsPdf *pdftmp = simPdf->getPdf(channelCat.getCurrentLabel());
      assert(pdftmp != 0);

      if (printLevel > 1) {
         std::cout << "on type " << channelCat.getCurrentLabel() << " " << channelCat.getCurrentIndex() << std::endl;
      }

      RooDataSet *dataSinglePdf = generateAsimovSinglePdf(*pdftmp, observables, *weightVar, &channelCat, printLevel);
      if (!dataSinglePdf) {
         oocoutE((TObject *)0, Generation)
            << "Error generating an Asimov data set for pdf " << pdftmp->GetName() << std::endl;
         return nullptr;
      }

      if (asimovDataMap.count(std::string(channelCat.getCurrentLabel())) != 0) {
         oocoutE((TObject *)0, Generation)
            << "AsymptoticCalculator::GenerateAsimovData(): The PDF for " << channelCat.getCurrentLabel()
            << " was already defined. It will be overridden. The faulty category definitions follow:" << std::endl;
         channelCat.Print("V");
      }

      asimovDataMap[std::string(channelCat.getCurrentLabel())] = dataSinglePdf;

      if (printLevel > 1) {
         std::cout << "channel: " << channelCat.getCurrentLabel() << ", data: ";
         dataSinglePdf->Print();
         std::cout << std::endl;
      }
   }

   RooArgSet obsAndWeight(observables);
   obsAndWeight.add(*weightVar);

   auto asimovData = std::make_unique<RooDataSet>("asimovDataFullModel", "asimovDataFullModel",
                                                  RooArgSet(obsAndWeight, channelCat), RooFit::Index(channelCat),
                                                  RooFit::Import(asimovDataMap), RooFit::WeightVar(*weightVar));

   for (auto &element : asimovDataMap) {
      delete element.second;
   }

   return asimovData;
}
