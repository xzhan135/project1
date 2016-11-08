#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"

double timeScale = 1./60.;

void ScaleHistogram(TH1F* inputHisto, TH1F* outputHisto, double scale){
	outputHisto->Scale(0.);
	int nbin = inputHisto->GetNbinsX();
	for (int i = 0; i < nbin; i++){
		double scaledBinContent = inputHisto->GetBinContent(i)*scale;
		outputHisto->SetBinContent(i, scaledBinContent);
	}
	std::cout<<scale<<std::endl;
}

int main(){

//	TFile* f01 = new TFile("/home/xzhan135/data/P50Data/MC/Cf252-P50D-MC/Cf252-pw-2-1_Det_E.root");
//  f01->ls();

//	TH1F* hnew01 = new TH1F("hnew01", "MC Energy Deposition", 250, 0., 5.);
//	hnew01 = (TH1F*)f01->Get("E_dep_0");

	TFile* f02 = new TFile("/home/xzhan135/data/P50Data/MC/Cf252-P50D-MC/Cf252-pw-2-0_DetSim-Processed.root");
	TH1F* hnew02 = new TH1F("hnew02", "MC Energy Reconstructed", 250, 0., 5.);
	hnew02 = (TH1F*)f02->Get("hELi1");
	hnew02->Sumw2();

	TFile* f03 = new TFile("/home/xzhan135/data/P50Data/Data/Cf252-P50D-Data-S028-F016-Processed.root");
	TH1F* hnew03 = new TH1F("hnew03", "Data Energy Reconstructed", 250, 0., 5.);
	hnew03 = (TH1F*)f03->Get("hELi1");
	hnew03->Sumw2();	
	
	TFile* f04 = new TFile("/home/xzhan135/data/P50Data/Data/BKGD-P50D-Data-S029-F000-Processed.root");
	TH1F* hnew04 = new TH1F("hnew04", "Background Data", 250, 0., 5.);
	hnew04 = (TH1F*)f04->Get("hELi1");
	hnew04->Sumw2();	
	//hnew04->Scale(timeScale);

	TH1F* hclone04 = (TH1F*)hnew04->Clone("BKGD");
	
	ScaleHistogram(hnew04, hclone04, timeScale);
	
	//hnew03->Add(hclone04, -1);

  Double_t norm = 1;
//  hnew01->Scale(norm/hnew01->Integral(), "width");
  hnew02->Scale(norm/hnew02->Integral(), "width");
  hnew03->Scale(norm/hnew03->Integral(), "width");
	hclone04->Scale(norm/hclone04->Integral(), "width");

	TCanvas *c1 = new TCanvas("c1","Compare Detector Response",20,10,700,400);
  gStyle->SetOptStat(0);
  c1->cd();

  hnew02->SetTitle("P50D Cf252 at Middle, seg - 1");
	hnew02->GetXaxis()->SetRangeUser(0, 10.);
	hnew02->GetYaxis()->SetRangeUser(0, 5.);
  hnew02->SetLineWidth(2);
  hnew02->SetLineColor(kBlue);
  hnew02->Draw("SAME");
	hnew03->SetLineWidth(2);
  hnew03->SetLineColor(kGreen);
  hnew03->Draw("SAME");
	
  TLegend* l1 = new TLegend(0.75, 0.75, 0.95, 0.95);
  l1->SetFillColor(kWhite);
  l1->SetTextFont(132);
//  l1->AddEntry(hnew01, "MC E Deposition", "l");
  l1->AddEntry(hnew02, "MC E Reconstructed", "l");
  l1->AddEntry(hnew03, "Data E Reconstructedd", "l");
  l1->Draw();

  c1->SaveAs("CompareDetResponseNCaptLi-2-1.pdf");
	c1->SaveAs("CompareDetResponseNCaptLi-2-1.png");

  TFile* gOutputFile = new TFile("CompareDetResponseNCaptLi-2-1.root", "RECREATE");
  gOutputFile->cd();
	TH1F* hclone02 = (TH1F*)hnew02->Clone("MC_E_Spect");
	TH1F* hclone03 = (TH1F*)hnew03->Clone("Data_E_Spect");
  hclone02->Write();
  hclone03->Write();
	hclone04->Write();
  gOutputFile->Close();

  return 0;
}
