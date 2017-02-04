
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
double EMax = 10;
double ELiMin = 0.4;
double ELiMax = 0.7;
double ENoise = 0.05;
double YMax = 450.;
double histoWidth = 0.01;
int nBins = (EMax-EMin)/histoWidth;

int main(int argc, char** argv){
	if (argc < 2) std::cout<< "The input file is missing. \n";

	std::string fileList = argv[1];
	std::string InputFile;
  std::cout<<"Reading the data file: "<< fileList <<".."<<std::endl;

	std::string OutputFile = (fileList);
	std::size_t listPos = fileList.find(".list");
	OutputFile.replace(listPos, 5, "-Processed.root");	
	std::cout<<"check"<<std::endl;

	std::cout<< "Finished loading input file.. \n";

	TFile* gOutputFile = TFile::Open(OutputFile.c_str(), "RECREATE");

  // Make a PSD vs Energy histogram
  TH2F* hPSDvsE0 = new TH2F("hPSDvsE0", ";Energy(MeV);PSD", nBins, EMin, EMax, nBins, PSDMin, PSDMax);
  TH2F* hPSDvsE1 = new TH2F("hPSDvsE1", ";Energy(MeV);PSD", nBins, EMin, EMax, nBins, PSDMin, PSDMax);

  // Make a PSD histogram for the nLi peak
  TH1F* hPSDnLi0 = new TH1F("hPSDnLi0", "PSD at the nLi Peak; PSD; Events", nBins, PSDMin, PSDMax);
  TH1F* hPSDnLi1 = new TH1F("hPSDnLi1", "PSD at the nLi Peak; PSD; Events", nBins, PSDMin, PSDMax);

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
	TH1F* hTime0 = new TH1F("hTime0", "Time distribution of neutron recoil", 200, -200, 200);
	TH1F* hTime1 = new TH1F("hTime1", "Time distribution of neutron recoil", 200, -200, 200);

 // TH1F* hENCapTot = new TH1F("hENCapTot", "Other N-Cap Energy: Total; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap1 = new TH1F("hENCap1", "Other N-Cap Energy: Top Segment; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap0 = new TH1F("hENCap0", "Other N-Cap Energy: Bottom Segment; Energy (MeV); Events", nBins, EMin, EMax);

	int limitEntries = MCpptree->GetEntries_PP_Tree();
	Long64_t lastEvent = 0;
	double totEnergy;
	double totRunTime = MCpptree->gPhysPulse_RunTime;

	std::cout << "Reading " << limitEntries <<" pulses.. \n";
	std::cout << "Total run time is "<< totRunTime << "s.. \n";

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
		
		if (Segment == 0 || Segment == 1);
		totEnergy += Energy;
		


		if (Segment == 0) hPSDvsE0->Fill(Energy, PSD);
		if (Segment == 1) hPSDvsE1->Fill(Energy, PSD);
		if (Energy < ELiMax && Energy > ELiMin && Segment == 0) hPSDnLi0->Fill(PSD);
		if (Energy < ELiMax && Energy > ELiMin && Segment == 1) hPSDnLi1->Fill(PSD);
	
    if (Energy < ENoise) continue;
  }
	int haha = 0; 
	for (int ev = 0; ev < limitEntries; ev++){
		
		MCpptree->GetEntry_PP_Tree(ev);

		Long64_t Event = MCpptree->GetPPEvent();
		float Energy = MCpptree->GetPPEnergy();
		float PSD = MCpptree->GetPPPSD();
		double Time = MCpptree->GetPPTime();
		int Segment = MCpptree->GetPPSegment();
		float YPos = MCpptree->GetPPY();

		if (lastEvent != Event){
			lastEvent = Event;
			if (totEnergy > 0.) hEGammaTot->Fill(totEnergy);
			totEnergy = 0.0;
		}
		
		if (Segment == 0 || Segment == 1);
		totEnergy += Energy;			

    if (Energy < ENoise) continue;

		hPSDnLi0->GetXaxis()->SetRangeUser(0.15, 0.25);
		double PSDGMax0 = hPSDnLi0->GetBinCenter(hPSDnLi0->GetMinimumBin());
		hPSDnLi0->GetXaxis()->SetRangeUser(PSDMin, PSDMax);
   	
		hPSDnLi1->GetXaxis()->SetRangeUser(0.15, 0.25);
		double PSDGMax1 = hPSDnLi1->GetBinCenter(hPSDnLi1->GetMinimumBin());
		hPSDnLi1->GetXaxis()->SetRangeUser(PSDMin, PSDMax);

    // Cut out neutron-like sources
    	if (Segment == 0 && PSD < PSDGMax0) hEGamma0->Fill(Energy);
    	if (Segment == 0 && PSD < PSDGMax0) hYGamma0->Fill(YPos);

    	if (Segment == 1 && PSD < PSDGMax1) hEGamma1->Fill(Energy);
    	if (Segment == 1 && PSD < PSDGMax1) hYGamma1->Fill(YPos);

// Fill the time distribution
		// Set PSD upper limit.		
    if (PSD < 0.35) {
			// Set energy threshold for recoil proton events.
			if (Energy < 12 && Energy > 1){
				// Set volume specifically to 0, and PSD to be nuclear recoil.
				if (Segment == 0 && PSD > PSDGMax0)	{		
						double StartTime = Time;																				//Find time of the nuclear recoil;
						double StartPos = YPos;																				//Find position of the nuclear recoil;
						if (std::abs(StartPos)>1750) continue;
						//std::cout << "Reading event No. " << Event << ", energy: " << Energy<< " MeV, PSD: " << PSD << ", Segment: " << Segment << ", Y: " << YPos << " \n";
						// Re-read through the file to catch neutron candidates.
						haha++;
						int counti0 = 0;
						for (int evi; evi < limitEntries; evi++){
							MCpptreei0->GetEntry_PP_Tree(evi);

							Long64_t Eventi0 = MCpptreei0->GetPPEvent();
							float Energyi0 = MCpptreei0->GetPPEnergy();
							float PSDi0 = MCpptreei0->GetPPPSD();
							double Timei0 = MCpptreei0->GetPPTime();
							int Segmenti0 = MCpptreei0->GetPPSegment();
							float YPosi0 = MCpptreei0->GetPPY();

							if (Eventi0 == Event) continue;
							//std::cout<<"check1"<<std::endl;
							if (Segmenti0 != Segment) continue;
							//std::cout<<"check2"<<std::endl;
							if (PSDi0 > PSDGMax0) continue;
							//std::cout<<"check3"<<std::endl;
							if (std::abs(YPosi0-YPos) > 66.66) continue;
							//std::cout<<"check4"<<std::endl;
							if (std::abs(Timei0-Time) > 200000000) continue;
							//std::cout<<"check5"<<std::endl;
							counti0 ++;
							
							hTime0->Fill((Timei0-Time)/1000000);		
							
						}
						std::cout<<counti0<<std::endl;
					}
				if (Segment == 1 && PSD > PSDGMax1) hELi1->Fill(Energy);
			}

		}
  }
	
	std::cout<<haha<<std::endl;
	hEGammaTot->Sumw2();
	hEGamma0->Sumw2();
	hEGamma1->Sumw2();
	hELi0->Sumw2();
	hELi1->Sumw2();
	
	hEGammaTot->Scale(1./totRunTime);
	hEGamma0->Scale(1./totRunTime);
	hEGamma1->Scale(1./totRunTime);
	hELi0->Scale(1./totRunTime);
	hELi1->Scale(1./totRunTime);
  
  hPSDvsE0->Write();
  hPSDvsE1->Write();
  hPSDnLi0->Write();
  hPSDnLi1->Write();
  hEGammaTot->Write();
  hEGamma0->Write();
  hEGamma1->Write();
  hYGamma0->Write();
  hYGamma1->Write();
	hELi0->Write();
	hELi1->Write();
	hTime0->Write();
	//hENCap0->Write();
	//hENCap1->Write();
  
  gOutputFile->Write();
  gOutputFile->Close();

	return 0;
}
