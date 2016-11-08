#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"

#include "SetupDataTree.hh"

double PSDMin = 0.;
double PSDMax = 0.4;
double EMin = 0.;
double EMax = 5;
double ELiMin = 0.4;
double ELiMax = 0.7;
double ENoise = 0.05;
double PSDGMax = 0.2;
double YMax = 450.;
int nBins = 250;

int main(int argc, char** argv){
	if (argc < 2) std::cout<< "The input file is missing. \n";

	std::string InputFile = argv[1];
  std::cout<<"Reading the data file: "<< InputFile <<".."<<std::endl;

	std::string OutputFile = (InputFile);
	std::size_t listPos = InputFile.find(".root");
	OutputFile.replace(listPos, 5, "-Processed.root");	

	PhysPulseTree* MCpptree = new PhysPulseTree(InputFile);

	std::cout<< "Finished loading input file.. \n";

	TFile* gOutputFile = TFile::Open(OutputFile.c_str(), "RECREATE");

  // Make a PSD vs Energy histogram
  TH2F* hPSDvsE = new TH2F("hPSDvsE", ";Energy(MeV);PSD", nBins, EMin, EMax, nBins, PSDMin, PSDMax);

  // Make a PSD histogram for the nLi peak
  TH1F* hPSDnLi = new TH1F("hPSDnLi", "PSD at the nLi Peak; PSD; Events", nBins, PSDMin, PSDMax);

  // Make an Energy histogram of gamma like events
  TH1F* hEGammaTot = new TH1F("hEGammaTot", "Gamma-like Energy: Total; Energy (MeV); Events", nBins, EMin, EMax);
  TH1F* hEGamma1 = new TH1F("hEGamma1", "Gamma-like Energy: Top Segment; Energy (MeV); Events", nBins, EMin, EMax);
  TH1F* hEGamma0 = new TH1F("hEGamma0", "Gamma-like Energy: Bottom Segment; Energy (MeV); Events", nBins, EMin, EMax);
  
  // Make an y position histogram of gamma like events
  TH1F* hYGamma1 = new TH1F("hYGamma1", "Gamma-like Y Position: Top Segment; Y Position (mm); Events", nBins,  -1.*YMax, YMax);
  TH1F* hYGamma0 = new TH1F("hYGamma0", "Gamma-like Y Position: Bottom Segment; Y Position (mm); Events", nBins,  -1.*YMax, YMax);
  
 // TH1F* hELiTot = new TH1F("hELiTot", "Li N-Cap Energy: Total; Energy (MeV); Events", nBins, EMin, EMax);
  TH1F* hELi1 = new TH1F("hELi1", "Li N-Cap Energy: Top Segment; Energy (MeV); Events", nBins, EMin, EMax);
  TH1F* hELi0 = new TH1F("hELi0", "Li N-Cap Energy: Bottom Segment; Energy (MeV); Events", nBins, EMin, EMax);

 // TH1F* hENCapTot = new TH1F("hENCapTot", "Other N-Cap Energy: Total; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap1 = new TH1F("hENCap1", "Other N-Cap Energy: Top Segment; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap0 = new TH1F("hENCap0", "Other N-Cap Energy: Bottom Segment; Energy (MeV); Events", nBins, EMin, EMax);

	int limitEntries = MCpptree->GetEntries_PP_Tree();
	Long64_t lastEvent = 0;
	double totEnergy;

	std::cout << "Reading " << limitEntries <<" pulses.. \n";

	for (int ev = 0; ev < limitEntries; ev++){
		
		MCpptree->GetEntry_PP_Tree(ev);

		Long64_t Event = MCpptree->GetPPEvent();
		float Energy = MCpptree->GetPPEnergy();
		float PSD = MCpptree->GetPPPSD();
		int Segment = MCpptree->GetPPSegment();
		float YPos = MCpptree->GetPPY();

		if (Energy > EMax && Energy < EMin) continue;
		if (PSD > PSDMax && PSD < PSDMin) continue;
		if (Segment != 0 && Segment != 1) continue;

		if (lastEvent != Event){
			lastEvent = Event;
			if (totEnergy > 0.) hEGammaTot->Fill(totEnergy);
			totEnergy = 0.0;
		}
		else{
			if (Segment == 0 || Segment == 1);
			totEnergy += Energy;
		}

		if (ev < 25) std::cout << "Reading event No. " << Event << ", energy: " << Energy<< " MeV, PSD: " << PSD << ", Segment: " << Segment << ", Y: " << YPos << " \n";

		hPSDvsE->Fill(Energy, PSD);
		if (Energy < ELiMax && Energy > ELiMin) hPSDnLi->Fill(PSD);
	
    if (Energy < ENoise) continue;

		
    
    // Cut out neutron-like sources
    if (PSD < PSDGMax) {
    	if (Segment == 0) hEGamma0->Fill(Energy);
    	if (Segment == 0) hYGamma0->Fill(YPos);

    	if (Segment == 1) hEGamma1->Fill(Energy);
    	if (Segment == 1) hYGamma1->Fill(YPos);
		}

    if (PSD > PSDGMax && PSD < 0.35) {
//			if (Energy < ELiMax && Energy > ELiMin) {
				if (Segment == 0) hELi0->Fill(Energy);
				if (Segment == 1) hELi1->Fill(Energy);
	//		}

    //  else {
    //    if (Segment == 0) hENCap0->Fill(Energy);
    //    if (Segment == 1) hENCap1->Fill(Energy);
    //  }
		}
  }
  
  hPSDvsE->Write();
  hPSDnLi->Write();
  hEGammaTot->Write();
  hEGamma0->Write();
  hEGamma1->Write();
  hYGamma0->Write();
  hYGamma1->Write();
	hELi0->Write();
	hELi1->Write();
	//hENCap0->Write();
	//hENCap1->Write();
  
  gOutputFile->Write();
  gOutputFile->Close();

	return 0;
}
