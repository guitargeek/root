// @(#)root/roostats:$Id$
// Author: Kyle Cranmer   28/07/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "Rtypes.h"

#if !defined(R__ALPHA) && !defined(R__SOLARIS) && !defined(R__ACC) && !defined(R__FBSD)
NamespaceImp(RooStats)
#endif

#include "TTree.h"
#include "TBranch.h"

#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooUniform.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"

using std::endl, std::vector;


namespace {
     template<class listT, class stringT> void getParameterNames(const listT* l,std::vector<stringT>& names){
       // extract the parameter names from a list
       if(!l) return;
       for (auto const *obj : *l) {
         names.push_back(obj->GetName());
       }
     }
     void getArgs(RooWorkspace* ws, const std::vector<TString> names, RooArgSet& args){
       for(const auto& p:names){
         RooRealVar* v = ws->var(p.Data());
         if(v){
           args.add(*v);
         }
       }
     }
   }

namespace RooStats {

   RooStatsConfig& GetGlobalRooStatsConfig() {
      static RooStatsConfig theConfig;
      return theConfig;
   }

   double AsimovSignificance(double s, double b, double sigma_b ) {
   // Asimov significance
   // formula [10] and [20] from  https://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf
      // case we have a sigma_b
      double sb2 = sigma_b*sigma_b;
      // formula below has a large error when sigma_b becomes zero
      // better to use the approximation for sigma_b=0 for very small values
      double r = sb2/b;
      if (r > 1.E-12) {
         double bpsb2 = b + sb2;
         double b2 = b*b;
         double spb = s+b;
         double za2 = 2.*( (spb)* std::log( ( spb)*(bpsb2)/(b2+ spb*sb2) ) -
                           (b2/sb2) * std::log(1. + ( sb2 * s)/(b * bpsb2) ) );
         return sqrt(za2);

      }
      // case when the background (b) is known
      double za2 = 2.*( (s+b) * std::log(1. + s/b) -s );
      return std::sqrt(za2);
   }

   RooAbsPdf * MakeNuisancePdf(const RooStats::ModelConfig &model, const char *name) {
      // The constraint terms are no part of the pdf of the model: the
      // ModelConfig keeps them as a separate set.
      const RooArgSet *constraints = model.GetExternalConstraints();
      if (!constraints || constraints->empty()) {
         oocoutW(nullptr, Eval) << "RooStatsUtils::MakeNuisancePdf - the input model has no constraint terms" << std::endl;
         return nullptr;
      }
      return new RooProdPdf(name, "", *constraints);
   }

   // Helper class for GetAsTTree
   class BranchStore {
      public:
         std::map<TString, double> fVarVals;
         double fInval;
         TTree *fTree = nullptr;

         BranchStore(const vector<TString> &params = vector<TString>(), double _inval = -999.) : fInval(_inval)
         {
         for (unsigned int i = 0; i < params.size(); i++)
            fVarVals[params[i]] = _inval;
         }

         ~BranchStore() {
            if (fTree) {
               for(std::map<TString, double>::iterator it = fVarVals.begin();it!=fVarVals.end();++it) {
                  TBranch *br = fTree->GetBranch( it->first );
                  if (br) br->ResetAddress();
               }
            }
         }

         void AssignToTTree(TTree &myTree) {
            fTree = &myTree;
            for(std::map<TString, double>::iterator it = fVarVals.begin();it!=fVarVals.end();++it) {
               const TString& name = it->first;
               myTree.Branch( name, &fVarVals[name], TString::Format("%s/D", name.Data()));
            }
         }
         void ResetValues() {
            for(std::map<TString, double>::iterator it = fVarVals.begin();it!=fVarVals.end();++it) {
               const TString& name = it->first;
               fVarVals[name] = fInval;
            }
         }
   };

   BranchStore* CreateBranchStore(const RooDataSet& data) {
      if (data.numEntries() == 0) {
         return new BranchStore;
      }
      vector <TString> V;
      for (auto *rvar : dynamic_range_cast<RooRealVar *>(* data.get(0))) {
         if (rvar == nullptr)
            continue;
         V.push_back(rvar->GetName());
         if (rvar->hasAsymError()) {
            V.push_back(TString::Format("%s_errlo", rvar->GetName()));
            V.push_back(TString::Format("%s_errhi", rvar->GetName()));
         }
         else if (rvar->hasError()) {
            V.push_back(TString::Format("%s_err", rvar->GetName()));
         }
      }
      return new BranchStore(V);
   }

   void FillTree(TTree &myTree, const RooDataSet &data) {
      BranchStore *bs = CreateBranchStore(data);
      bs->AssignToTTree(myTree);

      for(int entry = 0;entry<data.numEntries();entry++) {
         bs->ResetValues(); 
         for (auto const *rvar : dynamic_range_cast<RooRealVar *>(*data.get(entry))) {
            if (rvar == nullptr)
               continue;
            bs->fVarVals[rvar->GetName()] = rvar->getValV();
            if (rvar->hasAsymError()) {
               bs->fVarVals[TString::Format("%s_errlo", rvar->GetName())] = rvar->getAsymErrorLo();
               bs->fVarVals[TString::Format("%s_errhi", rvar->GetName())] = rvar->getAsymErrorHi();
            }
            else if (rvar->hasError()) {
               bs->fVarVals[TString::Format("%s_err", rvar->GetName())] = rvar->getError();
            }
         }
         myTree.Fill();
      }

      delete bs;
   }

   TTree * GetAsTTree(TString name, TString desc, const RooDataSet& data) {
      TTree* myTree = new TTree(name, desc);
      FillTree(*myTree, data);
      return myTree;
   }


   // useful function to print in one line the content of a set with their values
   void PrintListContent(const RooArgList & l, std::ostream & os ) {
      bool first = true;
      os << "( ";
      for (std::size_t i = 0; i< l.size(); ++i) {
         if (first) {
            first=false ;
         } else {
            os << ", " ;
         }
         l[i].printName(os);
         os << " = ";
         l[i].printValue(os);
      }
      os << ")\n";
   }

   // clone a workspace, copying all needed components and discarding all others
   // start off with the old workspace
   RooWorkspace* MakeReducedWorkspace(RooWorkspace *oldWS, const char *newName,
                                   bool copySnapshots, const char *mcname,
                                   const char *newmcname, bool copyData) {
      auto objects = oldWS->allGenericObjects();
      RooStats::ModelConfig *oldMC =
          dynamic_cast<RooStats::ModelConfig *>(oldWS->obj(mcname));
      for (auto it : objects) {
        if (!oldMC) {
          oldMC = dynamic_cast<RooStats::ModelConfig *>(it);
        }
      }
      if (!oldMC)
        throw std::runtime_error("unable to retrieve ModelConfig");

      RooAbsPdf *origPdf = oldMC->GetPdf();

      // start off with the old modelconfig
      std::vector<TString> poilist;
      std::vector<TString> nplist;
      std::vector<TString> obslist;
      std::vector<TString> globobslist;
      std::vector<TString> constraintlist;
      RooAbsPdf *pdf = nullptr;
      if (oldMC) {
        pdf = oldMC->GetPdf();
        ::getParameterNames(oldMC->GetParametersOfInterest(), poilist);
        ::getParameterNames(oldMC->GetNuisanceParameters(), nplist);
        ::getParameterNames(oldMC->GetObservables(), obslist);
        ::getParameterNames(oldMC->GetGlobalObservables(), globobslist);
        ::getParameterNames(oldMC->GetExternalConstraints(), constraintlist);
      }
      if (!pdf) {
        if (origPdf)
          pdf = origPdf;
      }
      if (!pdf) {
        return nullptr;
      }

      // create them anew
      RooWorkspace *newWS = new RooWorkspace(newName ? newName : oldWS->GetName());
      newWS->autoImportClassCode(true);
      RooStats::ModelConfig *newMC = new RooStats::ModelConfig(newmcname, newWS);

      // Copy snapshots
      if (copySnapshots) {
        for (auto* snap : static_range_cast<RooArgSet const*>(oldWS->getSnapshots())) {
          newWS->saveSnapshot(snap->GetName(), *snap, /*importValuesdf*/true);
        }
      }

      newWS->import(*pdf, RooFit::RecycleConflictNodes());
      RooAbsPdf *newPdf = newWS->pdf(pdf->GetName());
      newMC->SetPdf(*newPdf);

      // The constraint terms are no part of the pdf, so they are imported separately
      if (oldMC->GetExternalConstraints()) {
        for (RooAbsArg *constraint : *oldMC->GetExternalConstraints()) {
          newWS->import(*constraint, RooFit::RecycleConflictNodes());
        }
      }

      if (copyData) {
        for (auto d : oldWS->allData()) {
          newWS->import(*d);
        }
      }

      RooArgSet poiset;
      ::getArgs(newWS, poilist, poiset);
      RooArgSet npset;
      ::getArgs(newWS, nplist, npset);
      RooArgSet obsset;
      ::getArgs(newWS, obslist, obsset);
      RooArgSet globobsset;
      ::getArgs(newWS, globobslist, globobsset);

      newMC->SetParametersOfInterest(poiset);
      newMC->SetNuisanceParameters(npset);
      newMC->SetObservables(obsset);
      newMC->SetGlobalObservables(globobsset);
      if (!constraintlist.empty()) {
        RooArgSet constraintset;
        for (const auto &name : constraintlist) {
          if (RooAbsPdf *constraint = newWS->pdf(name.Data())) {
            constraintset.add(*constraint);
          }
        }
        newMC->SetExternalConstraints(constraintset);
      }
      newWS->import(*newMC);

      return newWS;
  }

}
