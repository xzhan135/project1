#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"

#include "SetupDetResponseTree.hh"

double PSDMin = 0.;
double PSDMax = 0.4;
double EMin = 0.;
double EMax = 2.0;
double ELiMin = 0.4;
double ELiMax = 0.7;
double ENoise = 0.2;
double PSDGMax = 0.2;
double YMax = 450.;
int nBins = 100;

int main(int argc, char** argv){
	if (argc < 2) std::cout<< "The input file is missing. \n";

	std::string InputFile = argv[1];
  std::cout<<"Reading the data file: "<< InputFile <<".."<<std::endl;

	std::string OutputFile = (InputFile);
	std::size_t listPos = InputFile.find(".root");
	OutputFile.replace(listPos, 5, "-Processed.root");	

	PhysPulseTree* MCpptree = new PhysPulseTree();
	
	MCpptree->SetupPhysPulseTree(InputFile);

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
  
	int limitEntries = MCpptree->GetEntries_PP_Tree();
	Long64_t lastEvent = 0;
	double totEnergy;

	std::cout << "Reading " << limitEntries <<" pulses.. \n";

	for (int ev = 0; ev < limitEntries; ev++){
		
		MCpptree->GetEntry_PP_Tree(ev);

		Long64_t Event = MCpptree->GetPPEvent();
		float Energy = MCpptree->GetPPEnergy();
		double PSD = MCpptree->GetPPPSD();
		int Segment = MCpptree->GetPPSegment();
		double YPos = MCpptree->GetPPY();

		if (Energy > EMax && Energy < EMin) continue;
		if (PSD > PSDMax && Energy < PSDMin) continue;
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
    if (PSD > PSDGMax) continue;
    
    if (Segment == 0) hEGamma0->Fill(Energy);
    if (Segment == 0) hYGamma0->Fill(YPos);

    if (Segment == 1) hEGamma1->Fill(Energy);
    if (Segment == 1) hYGamma1->Fill(YPos);
    
  }
  
  hPSDvsE->Write();
  hPSDnLi->Write();
  hEGammaTot->Write();
  hEGamma0->Write();
  hEGamma1->Write();
  hYGamma0->Write();
  hYGamma1->Write();
  
  gOutputFile->Write();
  gOutputFile->Close();

	return 0;
}
