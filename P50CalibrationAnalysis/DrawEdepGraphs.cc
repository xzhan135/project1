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

double timeScale = 1./6.;

void ScaleHistogram(TH1F* inputHisto, TH1F* outputHisto, double scale){
	outputHisto->Scale(0.);
	int nbin = inputHisto->GetNbinsX();
	for (int i = 0; i < nbin; i++){
		double scaledBinContent = inputHisto->GetBinContent(i)*scale;
		outputHisto->SetBinContent(i, scaledBinContent);
	}
	std::cout<<scale<<std::endl;
}

int main(int argc, char** argv){
  if (argc < 2) std::cout<< "The input file is missing. \n";

  std::string fileList = argv[1];
  std::string InputFile;
  std::cout<<"Reading the processed data and MC file: "<< fileList <<".."<<std::endl;

  std::string OutputFile = (fileList);

  std::size_t listPos = fileList.find(".list");
  OutputFile.replace(listPos, 5, ".root");
  std::cout<<"check"<<std::endl;

	std::ifstream inputList(fileList.c_str());
	std::string inputFile;
	std::vector<TFile*> vFiles;

	while(inputList >> inputFile){
		if (inputList.fail()) break;
		vFiles.push_back(new TFile(inputFile.c_str()));
	}
	
	vFiles.at(0)->ls();
	TH1F* hMC0 = new TH1F("hMC0", "MC Energy Deposition Seg 0", 500, 0., 5.);
	hMC0 = (TH1F*)vFiles.at(0)->Get("E_dep_0");
	TH1F* hMC1 = new TH1F("hMC1", "MC Energy Deposition Seg 1", 500, 0., 5.);
	hMC1 = (TH1F*)vFiles.at(0)->Get("E_dep_1");
	
	vFiles.at(1)->ls();
	TH1F* hData0 = new TH1F("hData0", "Data Energy Deposition Seg 0", 500, 0., 5.);
	hData0 = (TH1F*)vFiles.at(1)->Get("hEGamma0");
	TH1F* hLi0 = (TH1F*)vFiles.at(1)->Get("hELi0");
	TH1F* hData1 = new TH1F("hData1", "Data Energy Deposition Seg 1", 500, 0., 5.);
	hData1 = (TH1F*)vFiles.at(1)->Get("hEGamma1");
	TH1F* hLi1 = (TH1F*)vFiles.at(1)->Get("hELi1");

	vFiles.at(2)->ls();
	TH1F* hBKGD0 = new TH1F("hBKGD0", "BKGD Energy Deposition Seg 0", 500, 0., 5.);
	hBKGD0 = (TH1F*)vFiles.at(2)->Get("hEGamma0");
	TH1F* hLiBKGD0 = (TH1F*)vFiles.at(2)->Get("hELi0");
	TH1F* hBKGD1 = new TH1F("hBKGD1", "BKGD Energy Deposition Seg 1", 500, 0., 5.);
	hBKGD1 = (TH1F*)vFiles.at(2)->Get("hEGamma1");
	TH1F* hLiBKGD1 = (TH1F*)vFiles.at(2)->Get("hELi1");
	std::cout<<"Check 2"<<std::endl;

	hData0->Add(hBKGD0, -1);
	hData1->Add(hBKGD1, -1);
	hLi0->Add(hLiBKGD0, -1);
	hLi1->Add(hLiBKGD1, -1);
/*
	TCanvas *c1 = new TCanvas("c1","Compare Detector Response",20,10,700,400);
  gStyle->SetOptStat(0);
  c1->cd();

  hnew02->SetTitle("P50D Na22 at Middle, seg - 0");
	hnew02->GetXaxis()->SetRangeUser(0, 10.);
//	hnew02->GetYaxis()->SetRangeUser(0, 5.);
  hnew02->SetLineWidth(2);
  hnew02->SetLineColor(kBlue);
  hnew02->Draw("SAME");
	hnew03->SetLineWidth(2);
  hnew03->SetLineColor(kGreen);
  hnew03->Draw("SAME");
	hclone04->Draw("SAME");
	
  TLegend* l1 = new TLegend(0.75, 0.75, 0.95, 0.95);
  l1->SetFillColor(kWhite);
  l1->SetTextFont(132);
//  l1->AddEntry(hnew01, "MC E Deposition", "l");
  l1->AddEntry(hnew02, "MC E Reconstructed", "l");
  l1->AddEntry(hnew03, "Data E Reconstructedd", "l");
  l1->Draw();

  c1->SaveAs("CompareDetResponseNa22-2-0.pdf");
	c1->SaveAs("CompareDetResponseNa22-2-0.png");
*/

 	TFile* gOutputFile = new TFile(OutputFile.c_str(), "RECREATE");
  gOutputFile->cd();
	std::cout<<"check iii"<<std::endl;
	TH1F* hclone0 = (TH1F*)hMC0->Clone("MC_E_Spect_0");
	TH1F* hclone1 = (TH1F*)hData0->Clone("Data_E_Spect_0");
	std::cout<<"Check IV"<<std::endl;
	TH1F* hclone2 = (TH1F*)hBKGD0->Clone("BKGD_E_Spect_0");
	TH1F* hclone3 = (TH1F*)hMC1->Clone("MC_E_Spect_1");
	TH1F* hclone4 = (TH1F*)hData1->Clone("Data_E_Spect_1");
	TH1F* hclone5 = (TH1F*)hBKGD1->Clone("BKGD_E_Spect_1");
	TH1F* hclone6 = (TH1F*)hLi0->Clone("Data_Li_0");
	TH1F* hclone7 = (TH1F*)hLi1->Clone("Data_Li_1");
	hclone0->Write();
	hclone1->Write();
	hclone2->Write();
	hclone3->Write();
	hclone4->Write();
	hclone5->Write();
	hclone6->Write();
	hclone7->Write();
	gOutputFile->Write();
  gOutputFile->Close();

  return 0;
}
