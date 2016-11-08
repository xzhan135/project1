#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"

#include "DetResponse.hh"

DetResponse::DetResponse(){
  E_prim = new TH1F("E_prim", "Energy Spectrum of Primary Particles", 250, 0., 5.0);
  E_dep_tot = new TH1F("E_dep_tot", "Detected Energy Specturm", 250, 0., 5.);
	E_dep_0 = new TH1F("E_dep_0", "Energy Deposition, Seg 0", 250, 0.,5.);
	E_dep_1 = new TH1F("E_dep_1", "Energy Deposition, Seg 1", 250, 0.,5.);
}

void DetResponse::FillPrimEvent(ParticleVertex* ep, Long64_t ev){
  if (ep->PID < 65536){
    if (ev==testEvt) {
      energy += ep->E;
    }
    else {
      E_prim->Fill(energy);
      energy= ep->E;
      testEvt=ev;
    }
  }
}

void DetResponse::FillDetEvent(IoniCluster* ep, Long64_t ev){
  if (ev==testEvt) {
    energy += ep->Eq;
  }
  else {
		if (ep->vol == 1) E_dep_1->Fill(energy);
		if (ep->vol == 0) E_dep_0->Fill(energy);
    if (ep->vol == 0 || ep->vol == 1) E_dep_tot->Fill(energy);
    energy = ep->Eq;
    testEvt=ev;
  }
}

void DetResponse::DrawHisto(){
  TH1F* hnew = (TH1F*)E_prim->Clone("hnew");
  gStyle->SetOptStat(0);
  hnew->SetLineColor(3);
  hnew->GetXaxis()->SetTitle("Energy (MeV)");
  hnew->GetYaxis()->SetTitle("Entries");
  hnew->Draw("e1");

  TLegend* l1 = new TLegend(0.65, 0.15, 0.85, 0.35);
  l1->SetFillColor(kWhite);
  l1->SetTextFont(132);
  l1->AddEntry(hnew, "Primary Energy", "l");
  l1->Draw();

  TH1F* hnew1 = (TH1F*)E_dep_tot->Clone("hnew1");
  gStyle->SetOptStat(0);
  hnew1->SetLineColor(2);
  hnew1->GetXaxis()->SetTitle("Energy (MeV)");
  hnew1->GetYaxis()->SetTitle("Entries");
  hnew1->Draw("e1");

  TLegend* l2 = new TLegend(0.65, 0.15, 0.85, 0.35);
  l2->SetFillColor(kWhite);
  l2->SetTextFont(132);
  l2->AddEntry(hnew, "Primary Energy", "l");
  l2->Draw();
}

