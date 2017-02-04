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
#include "TPad.h"
#include "TLine.h"
#include "TPaveText.h"

double scaleLow = 0.9;
double scaleHigh = 1.1;
double scaleRes = 0.001;
double nbin_chi2 = (scaleHigh-scaleLow)/scaleRes;
double resLow = 0.02;
double resHigh = 0.1;
double resWidth = 0.001;
double nbin_res = (resHigh-resLow)/resWidth;

// Scale a input histogram to a target histogram with a user input scale..
void ScaleHistogram(TH1F* inputHist, TH1F* toyHist, double scale) {
  toyHist->Scale(0);
  
  // Create a temporary histogram.
  TH1F* tempHist = (TH1F*)inputHist->Clone("tempHist");
  
  double nBins = inputHist->GetNbinsX();
  
  for (int i = 0; i < nBins; i++) {
    // Clear temporary histogram
    tempHist->Scale(0);
    
    double content = inputHist->GetBinContent(i+1);
    double center = inputHist->GetBinCenter(i+1);
    double width = inputHist->GetBinWidth(i+1);
    double lowEdge = inputHist->GetBinLowEdge(i+1);
    double highEdge = lowEdge + width;
    
    // Find the bin that corresponds to the low edge scaled
    double newLowEdge = lowEdge*scale;
    double newLoBin = tempHist->FindBin(newLowEdge);
    
    // Find the bin that corresponds to the high edge scaled
    double newHighEdge = highEdge*scale;
    double newHiBin = tempHist->FindBin(newHighEdge);
    
    // Find the size of the new width for the content
    double newWidth = width * scale;
    
    // Loop through narrow bin range
		if (scale <= 1){
    	for (int j = newLoBin; j < newHiBin+1; j++) {
      	// Find Bin Edges
      	double low = inputHist->GetBinLowEdge(j);
      	double high = low+width;
      
      	double tempWidth = 0;
      
				if (low < newLowEdge && high > newHighEdge) tempWidth = newWidth;
    		else if (low < newLowEdge) tempWidth = high-newLowEdge;
      	else if (high > newHighEdge) tempWidth = newHighEdge - low;
      	else tempWidth = newWidth;
      
      	// Shift the histogram by a scaling factor
      	tempHist->SetBinContent(j, content * tempWidth / newWidth);
    	}
		}
	
		else if (scale > 1){
    	for (int j = newLoBin; j < newHiBin+1; j++) {
      	// Find Bin Edges
      	double low = inputHist->GetBinLowEdge(j);
      	double high = low+width;
      
      	double tempWidth = 0;
      
    		if (low < newLowEdge) tempWidth = high-newLowEdge;
      	else if (high > newHighEdge) tempWidth = newHighEdge - low;
      	else tempWidth = width;
      
      	// Shift the histogram by a scaling factor
      	tempHist->SetBinContent(j, content * tempWidth / newWidth);
    	}
		}
  	toyHist->Add(tempHist);
  }
}

// Compare the chi-square value between an input histogram to a reference histogram.
// This function allow user to set energy range. It will declare an in-function histogram
// for each of the two histogram with respect to the range, then compare the two weighted
// duplicated histogram by TH1F::Chi2Test(). The function returns chi-square value.
double Chi2Value(TH1F* inputHisto, TH1F* refHisto, double startBinContent, double endBinContent){
	int startBinNo = refHisto->FindBin(startBinContent);
	int endBinNo = refHisto->FindBin(endBinContent);
	int nBin = endBinNo-startBinNo;
	TH1F* hnew1 = (TH1F*)inputHisto->Clone("hnew1");
	TH1F* hnew2 = (TH1F*)refHisto->Clone("hnew2");
	double chi2;
	//for (int i = 1; i<=nBin; i++){
	//	hnew1->SetBinContent(i, inputHisto->GetBinContent(i+startBinNo));
	//	hnew2->SetBinContent(i, refHisto->GetBinContent(i+startBinNo));
	//}
	hnew1->GetXaxis()->SetRangeUser(startBinContent,endBinContent); 
	hnew2->GetXaxis()->SetRangeUser(startBinContent,endBinContent); 
	chi2 = hnew1->Chi2Test(hnew2, "Chi2/NDF WW");
	delete hnew1;
	delete hnew2;
	if (std::isinf(chi2)) chi2 = 10000;
	return chi2;
}

/*double Chi2Value(TH1F* inputHisto, TH1F* refHisto, double startBinContent, double endBinContent){
	int startBinNo = refHisto->FindBin(startBinContent);
	int endBinNo = refHisto->FindBin(endBinContent);
	double chi2;
	//std::cout<<"start at bin "<<startBinNo<<"; end at bin "<<endBinNo<<std::endl;
	for (int i = startBinNo; i<endBinNo; i++){
		double variance = pow((inputHisto->GetBinContent(i) - refHisto->GetBinContent(i)), 2);
		double expect = refHisto->GetBinContent(i);
		//std::cout<<inputHisto->GetBinContent(i) - refHisto->GetBinContent(i)<<std::endl;
		chi2 += (variance)/expect;
	}
	if (std::isinf(chi2)) chi2 = 10000;
	return chi2;
}*/

void GausSmear(TH1F* inputHist, double res){
  if(!inputHist) return;
  TH1F* hnew = (TH1F*)inputHist->Clone("hnew");
  inputHist->Reset();
  unsigned int nBins = inputHist->GetNbinsX();

  for (unsigned int i = 0; i < nBins; i++){
    double x = hnew->GetBinCenter(i+1);
		double sigma = res*TMath::Sqrt(x);
    double bandMin = hnew->FindBin(x - 5*sigma);
    double bandMax = hnew->FindBin(x + 5*sigma);
    double binContentSmeared = 0;
    
    for (unsigned int j = bandMin; j < bandMax; j++){
      double mean = hnew->GetBinCenter(j+1);
			double sigmaj = res*TMath::Sqrt(mean);
      binContentSmeared += hnew->GetBinContent(j+1)*TMath::Gaus(x,mean, sigmaj,false);
    }
    binContentSmeared=binContentSmeared*hnew->GetBinWidth(i+1)/TMath::Sqrt(2.*TMath::Pi())/sigma;
    inputHist->SetBinContent(i+1,binContentSmeared);
  }
  delete hnew;
}
	
void NormHistogram(TH1F* inputHist, TH1F* refHist, double startE, double endE){
	inputHist->Sumw2();
	//refHist->Sumw2();
	double startRange = inputHist->GetBinCenter(1);
	double endRange = inputHist->GetBinCenter(inputHist->GetNbinsX());
	inputHist->GetXaxis()->SetRange(inputHist->FindBin(startE), inputHist->FindBin(endE));
	double inputMax = inputHist->GetMaximum();
	refHist->GetXaxis()->SetRange(refHist->FindBin(startE),refHist->FindBin(endE));
	double refMax = refHist->GetMaximum();
	inputHist->Scale(refMax/inputMax);
	inputHist->GetXaxis()->SetRangeUser(startRange, endRange);
	refHist->GetXaxis()->SetRangeUser(startRange, endRange);
//	std::cout<<startRange<<"; "<< endRange<<std::endl;
}

void FillProjectionX(TH2* input2DHist, TH1* outputHist, int fixedYBin){
	int nBins = input2DHist->GetNbinsX();
	if (outputHist->GetNbinsX() != nBins) std::cout<< "h1 contains different X Channels! \n";
	for (int i = 1; i<= nBins; i++){
		outputHist->SetBinContent(i, input2DHist->GetBinContent(input2DHist->GetBin(i, fixedYBin)));
	}
}

void FillProjectionY(TH2* input2DHist, TH1* outputHist, int fixedXBin){
	int nBins = input2DHist->GetNbinsY();
	if (outputHist->GetNbinsX() != nBins) std::cout<< "h1 contains different X Channels! \n";
	for (int i = 1; i<= nBins; i++){
		outputHist->SetBinContent(i, input2DHist->GetBinContent(input2DHist->GetBin(fixedXBin, i)));
	}
}

int main(int argc, char** argv){
	
	std::string inputFile = argv[1];
	double startE = std::stod(argv[2]);
	double endE = std::stod(argv[3]);
		
	std::cout<<"Loading input file: "<< inputFile << ".. \n"; 
	
	std::string outputFile(inputFile);
	std::string matchedPDF_0(inputFile);
	std::string matchedPNG_0(inputFile);
	std::string chi2PDF_0(inputFile);
	std::string chi2PNG_0(inputFile);
	std::string matchedPDF_1(inputFile);
	std::string matchedPNG_1(inputFile);
	std::string chi2PDF_1(inputFile);
	std::string chi2PNG_1(inputFile);

	std::size_t listPos = inputFile.find(".root");

	if (endE < 1){
		outputFile.replace(listPos, 5, "-Matched-Low.root");
		matchedPDF_0.replace(listPos, 5, "-Matched-Low_0.pdf");
		matchedPNG_0.replace(listPos, 5, "-Matched-Low_0.png");
		chi2PDF_0.replace(listPos, 5, "-Chi2-Low_0.pdf");
		chi2PNG_0.replace(listPos, 5, "-Chi2-Low_0.png");
		matchedPDF_1.replace(listPos, 5, "-Matched-Low_1.pdf");
		matchedPNG_1.replace(listPos, 5, "-Matched-Low_1.png");
		chi2PDF_1.replace(listPos, 5, "-Chi2-Low_1.pdf");
		chi2PNG_1.replace(listPos, 5, "-Chi2-Low_1.png");
	}

	else if (endE >= 1){
		outputFile.replace(listPos, 5, "-Matched-Hi.root");
		matchedPDF_0.replace(listPos, 5, "-Matched-Hi_0.pdf");
		matchedPNG_0.replace(listPos, 5, "-Matched-Hi_0.png");
		chi2PDF_0.replace(listPos, 5, "-Chi2-Hi_0.pdf");
		chi2PNG_0.replace(listPos, 5, "-Chi2-Hi_0.png");
		matchedPDF_1.replace(listPos, 5, "-Matched-Hi_1.pdf");
		matchedPNG_1.replace(listPos, 5, "-Matched-Hi_1.png");
		chi2PDF_1.replace(listPos, 5, "-Chi2-Hi_1.pdf");
		chi2PNG_1.replace(listPos, 5, "-Chi2-Hi_1.png");
	}

	std::cout<< "Setting output file: " << outputFile << ".." <<std::endl;
	
	TFile* f01 = new TFile(inputFile.c_str(), "READ");
	f01->ls();

	TH1F*	hMC0 = (TH1F*)f01->Get("MC_E_Spect_0");
	TH1F* hData0 = (TH1F*)f01->Get("Data_E_Spect_0");
	
//	std::cout<< hMC0->GetNbinsX() << "; " << hData0->GetNbinsX() <<". \n";
	TH1F*	hMC1 = (TH1F*)f01->Get("MC_E_Spect_1");
	TH1F* hData1 = (TH1F*)f01->Get("Data_E_Spect_1");

	//TH1F* hSmeared0 = (TH1F*)hMC0->Clone("MC_Smeared_0");
	//TH1F* hSmeared1 = (TH1F*)hMC1->Clone("MC_Smeared_1");

	TH1F* hChi0 = new TH1F("hChi0", "#Chi^{2} vs Scale ", nbin_chi2, scaleLow, scaleHigh);
	TH1F* hChi1 = new TH1F("hChi1", "#Chi^{2} vs Scale ", nbin_chi2, scaleLow, scaleHigh);

	TH1F* hERes0 = new TH1F("hERes0", "#Chi^{2} vs E Resolution", nbin_res, resLow, resHigh); 
	std::cout<<"check "<<nbin_res<<std::endl;
	TH1F* hERes1 = new TH1F("hERes1", "#Chi^{2} vs E Resolution", nbin_res, resLow, resHigh); 
	
	TH2F* hChi2D_0 = new TH2F("hChi2D_0", "#Chi^{2} vs E Resolution vs Scale",  nbin_res, resLow, resHigh, nbin_chi2, scaleLow, scaleHigh);
	TH2F* hChi2D_1 = new TH2F("hChi2D_1", "#Chi^{2} vs E Resolution vs Scale",  nbin_res, resLow, resHigh, nbin_chi2, scaleLow, scaleHigh);

	//GausSmear(hSmeared0, 0.04);
	//GausSmear(hSmeared1, 0.04);

	TH1F* hMatched0 = (TH1F*)hMC0->Clone("MC_Matched_0");
	TH1F* hMatched1 = (TH1F*)hMC1->Clone("MC_Matched_1");

	std::cout<<"Finished loading input file .. \n";	

	for (int j = 1; j <= nbin_res; j++){
		double smearRes = resLow+j*resWidth;
//		std::cout<< "Res = " << smearRes<< std::endl;
		TH1F* hSmeared0 = (TH1F*)hMC0->Clone("MC_Smeared_0");
		TH1F* hSmeared1 = (TH1F*)hMC1->Clone("MC_Smeared_1");
		GausSmear(hSmeared0, smearRes);
		GausSmear(hSmeared1, smearRes);
		for (int i = 1; i <= nbin_chi2; i++){
			double matchScale = scaleLow+i*scaleRes;
//			std::cout<< "Scale = " << matchScale << std::endl;
			TH1F* hShift0 = (TH1F*)hSmeared0->Clone("hShift0");
			TH1F* hShift1 = (TH1F*)hSmeared1->Clone("hShift1");
			ScaleHistogram(hSmeared0, hShift0, matchScale);
			ScaleHistogram(hSmeared1, hShift1, matchScale);
			NormHistogram(hShift0, hData0, startE, endE);
			NormHistogram(hShift1, hData1, startE, endE);
			double chi2_0 = Chi2Value(hData0, hShift0, startE, endE);
			double chi2_1 = Chi2Value(hData1, hShift1, startE, endE);
			hChi2D_0->SetBinContent(j, i,chi2_0); 
			hChi2D_1->SetBinContent(j, i,chi2_1);
			//std::cout<<"Chi2 value = "<< chi2_0 <<".\n";
			delete hShift0, hShift1;
		}
		delete hSmeared0, hSmeared1;
	}
//	hChi0->Sumw2();
//	hChi1->Sumw2();
	
	// Fill the 1D projection of the 2D histogram.
	int chi2_minxBin0;
	int chi2_minyBin0;
	int chi2_zbin0;
	hChi2D_0->GetBinXYZ(hChi2D_0->GetMinimumBin(), chi2_minxBin0, chi2_minyBin0, chi2_zbin0); // Look for the x and y bin No. of the 2D histogram;
	std::cout << hChi2D_0->GetMinimumBin()<<"; "<<hChi2D_0->GetBin(chi2_minxBin0,chi2_minyBin0) <<"; "<<chi2_minxBin0 <<"; "<<chi2_minyBin0<<std::endl;
	// Fill the X projection, with fixed Y value;
	FillProjectionX(hChi2D_0, hERes0, chi2_minyBin0);	
	FillProjectionY(hChi2D_0, hChi0, chi2_minxBin0);
	
	double bestERes_0 = hERes0->GetBinCenter(hERes0->GetMinimumBin());
	double bestScale_0 = hChi0->GetBinCenter(hChi0->GetMinimumBin());
	double bestChi2_0 = hChi0->GetBinContent(hChi0->GetMinimumBin());
	double chi2Sigma_0 = bestChi2_0+1.;
	std::cout<< "Best chi2 value for segment 0 = "<<hChi2D_0->GetMinimum()<< ". \n";
	std::cout<< "Best chi2 value for segment 0 = "<<hChi0->GetBinContent(hChi0->GetMinimumBin())<<", coresponding scale = "<< bestScale_0 << ". \n";
	std::cout<< "Best chi2 value for segment 0 = "<<hERes0->GetBinContent(hERes0->GetMinimumBin())<<", coresponding E resolution = "<< bestERes_0 << ". \n";
	TH1F* hSmeared0 = (TH1F*)hMC0->Clone("MC_Smeared_0");
	GausSmear(hSmeared0, bestERes_0);
	ScaleHistogram(hSmeared0, hMatched0, bestScale_0);
	NormHistogram(hMatched0, hData0, startE, endE);
	NormHistogram(hSmeared0, hData0, startE, endE);

	//std::cout<<"Original chi-square value: " << hSmeared0->Chi2Test(hMatched0, "WW CHI2")<< ". \n";
	
	// Fill the 1D projection of the 2D histogram.
	int chi2_minxBin1;
	int chi2_minyBin1;
	int chi2_zbin1;
	hChi2D_1->GetBinXYZ(hChi2D_1->GetMinimumBin(), chi2_minxBin1, chi2_minyBin1, chi2_zbin1); // Look for the x and y bin No. of the 2D histogram;
	std::cout << hChi2D_1->GetMinimumBin()<<"; "<<hChi2D_1->GetBin(chi2_minxBin1,chi2_minyBin1) <<"; "<<chi2_minxBin1 <<"; "<<chi2_minyBin1<<std::endl;
	// Fill the X projection, with fixed Y value;
	FillProjectionX(hChi2D_1, hERes1, chi2_minyBin1);	
	FillProjectionY(hChi2D_1, hChi1, chi2_minxBin1);
	
	double bestERes_1 = hERes1->GetBinCenter(hERes1->GetMinimumBin());
	double bestScale_1 = hChi1->GetBinCenter(hChi1->GetMinimumBin());
	double bestChi2_1 = hChi1->GetBinContent(hChi1->GetMinimumBin());
	double chi2Sigma_1 = bestChi2_1+1.;
	std::cout<< "Best chi2 value for segment 1 = "<<hChi2D_1->GetMinimum()<< ". \n";
	std::cout<< "Best chi2 value for segment 1 = "<<hChi1->GetBinContent(hChi1->GetMinimumBin())<<", coresponding scale = "<< bestScale_1 << ". \n";
	std::cout<< "Best chi2 value for segment 1 = "<<hERes1->GetBinContent(hERes1->GetMinimumBin())<<", coresponding E resolution = "<< bestERes_1 << ". \n";
	TH1F* hSmeared1 = (TH1F*)hMC1->Clone("MC_Smeared_1");
	GausSmear(hSmeared1, bestERes_1);
	ScaleHistogram(hSmeared1, hMatched1, bestScale_1);
	NormHistogram(hMatched1, hData1, startE, endE);
	NormHistogram(hSmeared1, hData1, startE, endE);

	TCanvas* c01 = new TCanvas("c01", "Compare Shifted Histograms", 20, 10, 700, 400);
	gStyle->SetOptStat(0);
	c01->cd();

	//hMatched0->SetTitle(("P50D " +analysisMode+ " at rod -" +std::to_string(calibNo)+ ", seg - " + std::to_string(segNo) + ", in E range ("+ argv[4]+", "+ argv[5]+") MeV").c_str());
	hData0->SetLineColor(kBlue);
	hData0->Draw();
	//hSmeared0->SetLineColor(kRed);
	//hSmeared0->Draw("SAME");
	hMC0->SetLineColor(kGreen);
	NormHistogram(hMC0, hData0, startE, endE);
  hData0->GetXaxis()->SetRangeUser(0,5);
	hMC0->Draw("SAME");
	hMatched0->SetLineColor(kRed);
	hMatched0->Draw("SAME");

  TLegend* l1 = new TLegend(0.7, 0.6, 0.9, 0.9);
  l1->SetFillColor(kWhite);
  l1->SetTextFont(132);
  l1->AddEntry(hMatched0, "MC E Reconstructed", "l");
  l1->AddEntry(hData0, "Data E Reconstructedd", "l");
	l1->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale_0) ).c_str(), "");
	l1->AddEntry((TObject*)0, ("E resolution = "+std::to_string(bestERes_0) ).c_str(), "");
	l1->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2_0) ).c_str(), "");
  l1->Draw();

	c01->SaveAs(matchedPDF_0.c_str());
	c01->SaveAs(matchedPNG_0.c_str());

	TCanvas* c02 = new TCanvas("c02", "Chi 2 figure", 20, 10, 900, 600);
	gStyle->SetOptStat(0);
	c02->cd();

 	TPad *center_pad02 = new TPad("center_pad02", "center_pad",0.3,0.4,1,1);
  center_pad02->Draw();

  TPad* left_pad02 = new TPad("right_pad02", "right_pad",0.0,0.4,0.3,1.0);
  left_pad02->Draw();

  TPad* bottom_pad02 = new TPad("bottom_pad02", "bottom_pad",0.3,0.0,1.0,0.4);
  bottom_pad02->Draw();

  TPad* corner_pad02 = new TPad("corner_pad02", "corner_pad",0.0,0.0,0.28,0.36);
  corner_pad02->Draw();

	center_pad02->cd();
	center_pad02->SetLogz();
	center_pad02->SetLeftMargin(0);
	center_pad02->SetBottomMargin(0);
	hChi2D_0->GetXaxis()->SetLabelSize(0);
	hChi2D_0->GetYaxis()->SetLabelSize(0);
	hChi2D_0->GetZaxis()->SetTitle("#Chi^{2}/NDF ");
	hChi2D_0->GetZaxis()->SetTitleSize(0.06);
	hChi2D_0->GetZaxis()->SetTitleOffset(-0.8);
	hChi2D_0->Draw("colz");

	left_pad02->cd();
	left_pad02->SetLogx();
	left_pad02->SetRightMargin(0);
	left_pad02->SetLeftMargin(0.2);
	left_pad02->SetBottomMargin(0);
	hChi0->GetYaxis()->SetLabelSize(0.05);
	hChi0->GetYaxis()->SetLabelOffset(-0.1);
	hChi0->GetXaxis()->SetTitle("E Scale");
	hChi0->GetXaxis()->SetTitleSize(0.05);
	hChi0->GetXaxis()->SetTitleOffset(1.5);
	hChi0->GetXaxis()->SetLabelSize(0.05);
	hChi0->SetFillColor(kBlue-6);
	hChi0->Draw("hbar");

	bottom_pad02->cd();
	bottom_pad02->SetLogy();
	bottom_pad02->SetTopMargin(0);
	bottom_pad02->SetBottomMargin(0.2);
	bottom_pad02->SetLeftMargin(0);
	hERes0->GetYaxis()->SetLabelSize(0.055);
	hERes0->GetYaxis()->SetLabelOffset(-0.06);
	hERes0->GetXaxis()->SetTitle("E Resolution");
	hERes0->GetXaxis()->SetTitleSize(0.055);
	hERes0->GetXaxis()->SetLabelSize(0.055);
	hERes0->GetYaxis()->SetTitleOffset(0.15);
	hERes0->SetFillColor(kBlue-6);
	hERes0->Draw("bar");

	corner_pad02->cd();
	TPaveText *pt02 = new TPaveText(.1,.1,1,1);
	pt02->AddText(("Scale = "+std::to_string(bestScale_0) ).c_str());
	pt02->AddText(("E resolution = "+std::to_string(bestERes_0) ).c_str());
	pt02->AddText(("#Chi^{2} = "+std::to_string(bestChi2_0) ).c_str());
	pt02->Draw();
  /*TLegend* l2 = new TLegend(0.7, 0.8, 0.9, 0.9);
  l2->SetFillColor(kWhite);
  l2->SetTextFont(132);
	l2->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale_0) ).c_str(), "");
	l2->AddEntry((TObject*)0, ("E resolution = "+std::to_string(bestERes_0) ).c_str(), "");
	l2->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2_0) ).c_str(), "");
  l2->Draw();*/

	c02->SaveAs(chi2PDF_0.c_str());
	c02->SaveAs(chi2PNG_0.c_str());

	TCanvas* c11 = new TCanvas("c11", "Compare Shifted Histograms", 20, 10, 700, 400);
	gStyle->SetOptStat(0);
	c11->cd();

	hData1->SetLineColor(kBlue);
	hData1->Draw();
	//hSmeared0->SetLineColor(kRed);
	//hSmeared0->Draw("SAME");
	hMC1->SetLineColor(kGreen);
	NormHistogram(hMC1, hData1, startE, endE);
  hData1->GetXaxis()->SetRangeUser(0,5);
	hMC1->Draw("SAME");
	hMatched1->SetLineColor(kRed);
	hMatched1->Draw("SAME");

  TLegend* l11 = new TLegend(0.7, 0.6, 0.9, 0.9);
  l11->SetFillColor(kWhite);
  l11->SetTextFont(132);
  l11->AddEntry(hMatched0, "MC E Reconstructed", "l");
  l11->AddEntry(hData0, "Data E Reconstructedd", "l");
	l11->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale_1) ).c_str(), "");
	l11->AddEntry((TObject*)0, ("E resolution = "+std::to_string(bestERes_1) ).c_str(), "");
	l11->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2_1) ).c_str(), "");
  l11->Draw();

	c11->SaveAs(matchedPDF_1.c_str());
	c11->SaveAs(matchedPNG_1.c_str());

	TCanvas* c12 = new TCanvas("c12", "Chi 2 figure", 20, 10, 900, 600);
	gStyle->SetOptStat(0);
	c12->cd();

 	TPad *center_pad12 = new TPad("center_pad12", "center_pad",0.3,0.4,1,1);
  center_pad12->Draw();

  TPad* left_pad12 = new TPad("right_pad12", "right_pad",0.0,0.4,0.3,1.0);
  left_pad12->Draw();

  TPad* bottom_pad12 = new TPad("bottom_pad12", "bottom_pad",0.3,0.0,1.0,0.4);
  bottom_pad12->Draw();

  TPad* corner_pad12 = new TPad("corner_pad12", "corner_pad",0.0,0.0,0.28,0.36);
  corner_pad12->Draw();

	center_pad12->cd();
	center_pad12->SetLogz();
	center_pad12->SetLeftMargin(0);
	center_pad12->SetBottomMargin(0);
	hChi2D_1->GetXaxis()->SetLabelSize(0);
	hChi2D_1->GetYaxis()->SetLabelSize(0);
	hChi2D_1->GetZaxis()->SetTitle("#Chi^{2}/NDF ");
	hChi2D_1->GetZaxis()->SetTitleSize(0.06);
	hChi2D_1->GetZaxis()->SetTitleOffset(-0.8);
	hChi2D_1->Draw("colz");

	left_pad12->cd();
	left_pad12->SetLogx();
	left_pad12->SetRightMargin(0);
	left_pad12->SetLeftMargin(0.2);
	left_pad12->SetBottomMargin(0);
	hChi1->GetYaxis()->SetLabelSize(0.05);
	hChi1->GetYaxis()->SetLabelOffset(-0.1);
	hChi1->GetXaxis()->SetTitle("E Scale");
	hChi1->GetXaxis()->SetTitleSize(0.05);
	hChi1->GetXaxis()->SetTitleOffset(1.5);
	hChi1->GetXaxis()->SetLabelSize(0.05);
	hChi1->SetFillColor(kBlue-6);
	hChi1->Draw("hbar");

	bottom_pad12->cd();
	bottom_pad12->SetLogy();
	bottom_pad12->SetTopMargin(0);
	bottom_pad12->SetBottomMargin(0.2);
	bottom_pad12->SetLeftMargin(0);
	hERes1->GetYaxis()->SetLabelSize(0.055);
	hERes1->GetYaxis()->SetLabelOffset(-0.06);
	hERes1->GetXaxis()->SetTitle("E Resolution");
	hERes1->GetXaxis()->SetTitleSize(0.055);
	hERes1->GetXaxis()->SetLabelSize(0.055);
	hERes1->GetYaxis()->SetTitleOffset(0.15);
	hERes1->SetFillColor(kBlue-6);
	hERes1->Draw("bar");

	corner_pad12->cd();
	TPaveText *pt12 = new TPaveText(.1,.1,1,1);
	pt12->AddText(("Scale = "+std::to_string(bestScale_1) ).c_str());
	pt12->AddText(("E resolution = "+std::to_string(bestERes_1) ).c_str());
	pt12->AddText(("#Chi^{2} = "+std::to_string(bestChi2_1) ).c_str());
	pt12->Draw();

/*  TLegend* l12 = new TLegend(0.7, 0.7, 0.9, 0.9);
  l12->SetFillColor(kWhite);
  l12->SetTextFont(132);
	l12->AddEntry((TObject*)0, ("Scale = "+std::to_string(bestScale_1) ).c_str(), "");
	l12->AddEntry((TObject*)0, ("E resolution = "+std::to_string(bestERes_1) ).c_str(), "");
	l12->AddEntry((TObject*)0, ("#Chi^{2} = "+std::to_string(bestChi2_1) ).c_str(), "");
  l12->Draw();*/

	c12->SaveAs(chi2PDF_1.c_str());
	c12->SaveAs(chi2PNG_1.c_str());

	TFile* gOutputFile = new TFile(outputFile.c_str(), "RECREATE");
	gOutputFile->cd();
	hMC0->Write();
	hData0->Write();
	hMatched0->Write();
	hChi0->Write();
	hMC1->Write();
	hData1->Write();
	hMatched1->Write();
	hChi1->Write();
	gOutputFile->Write();
	gOutputFile->Close();
	return 0;
}
