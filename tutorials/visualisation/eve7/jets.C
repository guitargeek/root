/// \file
/// \ingroup tutorial_eve_7
///  This example display only points in web browser
///
/// \macro_code
///

#include "TRandom.h"
#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveJetCone.hxx>

namespace REX = ROOT::Experimental;

const Double_t kR_min = 240;
const Double_t kR_max = 250;
const Double_t kZ_d = 300;

void makeJets(int N_Jets, REX::REveElement *jetHolder)
{
   TRandom &r = *gRandom;

   for (int i = 0; i < N_Jets; i++) {
      auto jet = new REX::REveJetCone(Form("Jet_%d", i));
      jet->SetCylinder(2 * kR_max, 2 * kZ_d);
      jet->AddEllipticCone(r.Uniform(-0.5, 0.5), r.Uniform(0, TMath::TwoPi()), 0.1, 0.2);
      jet->SetFillColor(kPink - 8);
      jet->SetLineColor(kViolet - 7);

      jetHolder->AddElement(jet);
   }
}

void jets()
{
   auto eveMng = REX::REveManager::Create();

   REX::REveElement *jetHolder = new REX::REveElement("Jets");
   eveMng->GetEventScene()->AddElement(jetHolder);
   makeJets(7, jetHolder);

   eveMng->Show();
}
