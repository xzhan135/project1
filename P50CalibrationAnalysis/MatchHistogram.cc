#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

double scaleLow = 0.7;
double scaleHigh = 1.2;
double scaleRes = 0.001;
double nbin_chi2 = (scaleHigh-scaleLow)/scaleRes;

void ShiftHistogram(TH1F* inputHisto, TH1F* outputHisto, double scale){
	outputHisto->Scale(0);
	
	TH1F* hnew = (TH1F*)inputHisto->Clone("hnew");

	double nBins = inputHisto->GetNbinsX();

	for (int i = 0; i < nBins; i++){
		hnew->Scale(0);
		double content = inputHisto->GetBinContent(i+1);
		double center = inputHisto->GetBinCenter(i+1);
		double width = inputHisto->GetBinWidth(i+1);
		double lowEdge = inputHisto->GetBinLowEdge(i+1);
		double highEdge = lowEdge+width;
		
		double newLowEdge = lowEdge*scale;
		double newLowBin = hnew->FindBin(newLowEdge);
		
		double newHighEdge = highEdge*scale;
		double newHighBin = hnew->FindBin(newHighEdge);

		double newWidth = width*scale;

		for (int j = newLowBin; j < newHighBin+1; j++){
			double low = inputHisto->GetBinLowEdge(j);
			double high = low+width;
			
			double tempWidth = 0;
			if (low < newLowEdge) tempWidth = high - newLowEdge;
			else if (high > newHighEdge) tempWidth = newHighEdge - low;
			else tempWidth = width;

			hnew->SetBinContent(j, content * tempWidth/newWidth);
		}

		outputHisto->Add(hnew);
	}
}

double Chi2Value(TH1F* inputHisto, TH1F* refHisto, double startBinContent, double endBinContent){
	int startBinNo = refHisto->FindBin(startBinContent);
	int endBinNo = refHisto->FindBin(endBinContent);
	double chi2;
	for (int i = startBinNo; i<endBinNo; i++){
		double variance = pow((inputHisto->GetBinContent(i) - refHisto->GetBinContent(i)), 2);
		double expect = refHisto->GetBinContent(i);
		chi2 += (variance/expect);
	}
	if (std::isinf(chi2)) chi2 = 10000;
	return chi2/(endBinNo-startBinNo-1);
}

void AddRes(TH1F* inputHisto, TH1F* ResHist, double ResVal){
	ResHist->Scale(0.0);
	TH1F* hnew = (TH1F*)inputHisto->Clone("hnew");
	hnew->Scale(0.0);

	std::cout<<inputHisto->GetTitle()<<" have "<< hnew->GetNbinsX() <<" bins"<<std::endl;

	TF1 *ResolutionFunc = new TF1("gausres", "gaus(0)");
	std::cout<<"Resolution Function Created"<<std::endl;

	for(int i = 0; i<inputHisto->GetNbinsX(); i++){
		double ResSigma = ResVal*TMath::Sqrt(inputHisto->GetBinCenter(i+1));
		if (ResSigma > 0){
			double ResNorm = inputHisto->GetBinContent(i+1)*inputHisto->GetBinWidth(i+1);

			ResolutionFunc->SetParameters(ResNorm/(ResSigma*TMath::Sqrt(2*TMath::Pi())), inputHisto->GetBinCenter(i+1), ResSigma);

			int OverBins = 2*(ResSigma/inputHisto->GetBinWidth(i+1)+1);
			for (int j=i+1-OverBins; j<i+1+OverBins; j++){
				hnew->AddBinContent(j, ResolutionFunc->Eval(hnew->GetBinCenter(j)));
			}
		}
	}
	ResHist->Add(hnew);
	delete hnew;
}

int main(int argc, char** argv){
	
	std::string analysisMode = argv[1];
	int calibNo = atoi(argv[2]);
	int segNo = atoi(argv[3]);
	double startE = std::stod(argv[4]);
	double endE = std::stod(argv[5]);
		
	std::string inputFile = "CompareDetResponse" + analysisMode + "-" + std::to_string(calibNo) + "-" + std::to_string(segNo) + ".root";

	std::cout<<"Loading input file: "<< inputFile << ".. \n"; 
	
	std::string outputFile(inputFile);
	std::string matchedPDF(inputFile);
	std::string matchedPNG(inputFile);
	std::string chi2PDF(inputFile);
	std::string chi2PNG(inputFile);

	std::size_t listPos = inputFile.find(".root");

	if (endE < 1){
		outputFile.replace(listPos, 5, "-Matched-Low.root");
		matchedPDF.replace(listPos, 5, "-Matched-Low.pdf");
		matchedPNG.replace(listPos, 5, "-Matched-Low.png");
		chi2PDF.replace(listPos, 5, "-Chi2-Low.pdf");
		chi2PNG.replace(listPos, 5, "-Chi2-Low.png");
	}

	else if (endE >= 1){
		outputFile.replace(listPos, 5, "-Matched-Hi.root");
		matchedPDF.replace(listPos, 5, "-Matched-Hi.pdf");
		matchedPNG.replace(listPos, 5, "-Matched-Hi.png");
		chi2PDF.replace(listPos, 5, "-Chi2-Hi.pdf");
		chi2PNG.replace(listPos, 5, "-Chi2-Hi.png");
	}

	std::cout<< "Setting output file: " << outputFile << ".." <<std::endl;
	
	TFile* f01 = new TFile(inputFile.c_str(), "READ");
	f01->ls();
	TH1F* hnew01 = new TH1F("hnew01", "MC Energy Recontructed", 250, 0., 5.);
	TH1F* hnew02 = new TH1F("hnew02", "Data Energy Reconstructed", 250, 0., 5.);
	TH1F* hMatched = new TH1F("hMatchedHisto", "MC Matched with Data", 250, 0., 5.);
	TH1F* hChi2 = new TH1F("hChi2","Chi2 Value vs Scale Applied", nbin_chi2, scaleLow, scaleHigh);	
	TH1F* hSmeared = new TH1F("hSmeared", "MC Energy Smeared", 250, 0., 5.);

	hnew01 = (TH1F*)f01->Get("MC_E_Spect");
	hnew02 = (TH1F*)f01->Get("Data_E_Spect");
	
	std::cout<<"Finished loading input file .. \n";
	
	AddRes(hnew01, hSmeared, 0.035);

	for (int i = 0; i < nbin_chi2+1; i++){
		double matchScale = scaleLow+i*scaleRes;
		TH1F* hShift = new TH1F("hShift", "Shifted temporary histogram", 250, 0., 5.);
		ShiftHistogram(hSmeared, hShift, matchScale);
		double chi2 = Chi2Value(hnew02, hShift, startE, endE);
		hChi2->SetBinContent(i,chi2);
		//std::cout<<"Chi2 value = "<< chi2 <<".\n";
		delete hShift;
	}
	hChi2->Sumw2();
	
	double bestScale = hChi2->GetBinCenter(hChi2->GetMinimumBin());
	double bestChi2 = hChi2->GetBinContent(hChi2->GetMinimumBin());
	std::cout<< "Best chi2 value = "<<hChi2->GetBinContent(hChi2->GetMinimumBin())<<", coresponding scale = "<< bestScale << ". \n";
	ShiftHistogram(hSmeared, hMatched, bestScale);

	TCanvas* c01 = new TCanvas("c01", "Compare Shifted Histograms", 20, 10, 700, 400);
	gStyle->SetOptStat(0);
	c01->cd();

	hMatched->SetTitle(("P50D " +analysisMode+ " at rod -" +std::to_string(calibNo)+ ", seg - " + std::to_string(segNo) + ", in E range ("+ argv[4]+", "+ argv[5]+") MeV").c_str());
	hMatched->GetXaxis()->SetRangeUser(0,2.5);
	hMatched->SetLineColor(kRed);
	hMatched->Draw();
	hnew02->SetLineColor(kBlue);
	hnew02->Draw("SAME");

  TLegend* l1 = new TLegend(0.7, 0.6, 0.9, 0.9);
  l1->SetFillColor(kWhite);
  l1->SetTextFont(132);
  l1->AddEntry(hMatched, "MC E Reconstructed", "l");
  l1->AddEntry(hnew02, "Data E Reconstructedd", "l");
	l1->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale) ).c_str(), "");
	l1->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2) ).c_str(), "");
  l1->Draw();

	c01->SaveAs(matchedPDF.c_str());
	c01->SaveAs(matchedPNG.c_str());

	TCanvas* c02 = new TCanvas("c02", "Chi 2 figure", 20, 10, 700, 400);
	gStyle->SetOptStat(0);
	c02->cd();

	hChi2->SetTitle(("#Chi^{2} vs Scale ("+analysisMode+ " at rod -" +std::to_string(calibNo)+ ", seg - " + std::to_string(segNo)+", in E range ("+ argv[4]+", "+ argv[5]+") MeV)").c_str());
	//hChi2->GetYaxis()->SetRangeUser(0, 2);
	hChi2->Draw();

  TLegend* l2 = new TLegend(0.7, 0.8, 0.9, 0.9);
  l2->SetFillColor(kWhite);
  l2->SetTextFont(132);
	l2->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale) ).c_str(), "");
	l2->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2) ).c_str(), "");
  l2->Draw();

	c02->SaveAs(chi2PDF.c_str());
	c02->SaveAs(chi2PNG.c_str());

	TFile* gOutputFile = new TFile(outputFile.c_str(), "RECREATE");
	gOutputFile->cd();
	hnew01->Write();
	hnew02->Write();
	hMatched->Write();
	hChi2->Write();
	gOutputFile->Close();
	return 0;
}
