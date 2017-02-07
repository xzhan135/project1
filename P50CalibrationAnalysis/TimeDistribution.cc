
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TVector.h"

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
double PhysPulseRunTime = 0.0;

int main(int argc, char** argv){
	if (argc < 2) std::cout<< "The input file is missing. \n";

	std::string fileList = argv[1];
  std::ifstream inputList(fileList.c_str());
	std::string InputFile;
  std::cout<<"Reading the data file: "<< fileList <<".."<<std::endl;

	std::string OutputFile = (fileList);
	std::size_t listPos = fileList.find(".list");
	OutputFile.replace(listPos, 5, "-Processed.root");
	std::cout<<"check"<<std::endl;

	std::cout<< "Finished loading input file.. \n";


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
	TH1F* hTime0 = new TH1F("hTime0", "Time distribution of neutron recoil", 500, -100, 100);
	TH1F* hTime1 = new TH1F("hTime1", "Time distribution of neutron recoil", 200, -200, 200);

	//TH1F* hEB0 =new TH1F("hELi1", "Li N-Cap Energy: Top Segment; Energy (MeV); Events", nBins*2, EMin*2, EMax*2);
 // TH1F* hENCapTot = new TH1F("hENCapTot", "Other N-Cap Energy: Total; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap1 = new TH1F("hENCap1", "Other N-Cap Energy: Top Segment; Energy (MeV); Events", nBins, EMin, EMax);
 // TH1F* hENCap0 = new TH1F("hENCap0", "Other N-Cap Energy: Bottom Segment; Energy (MeV); Events", nBins, EMin, EMax);

	while (inputList >> InputFile) {
		if (inputList.fail()) break;
		if (InputFile[0] == '#') {
			continue;
		}
		if(InputFile.find("break") != std::string::npos) break;
		std::cout<<"Reading input file: "<< InputFile << ".." << std::endl;

		PhysPulseTree* MCpptree = new PhysPulseTree(InputFile);

		int limitEntries = MCpptree->GetEntries_PP_Tree();
		Long64_t lastEvent = 0;
		double totEnergy;

		std::cout << "Reading " << limitEntries <<" pulses.. \n";
		// Extract the PSD histograms add first.
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
	  }

		//Fill the time distribution histogram with nuclear recoil as start event and
		//electronic recoil as end event.
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
							if (std::abs(StartPos)>455) continue;
							std::cout << "Reading entry: "<<ev <<", event No. " << Event << ", energy: " << Energy<< " MeV, PSD: " << PSD << ", Segment: " << Segment << ", Y: " << YPos <<", Time: " << Time<< " \n";
							int formerEvt = ev-1;
							int latterEvt = ev+1;
							MCpptree->GetEntry_PP_Tree(formerEvt);
							double formerTime = MCpptree->GetPPTime();
							//std::cout<<"Former events to entry: " <<formerEvt <<", event: "<< MCpptree->GetPPEvent() <<std::endl;
							int formerCount = 0;
							while (formerEvt>=0 && (formerTime-StartTime)<0 && (formerTime-StartTime)>-200000) {
								double formerPos = MCpptree->GetPPY();
								double formerPSD = MCpptree->GetPPPSD();
								int formerSeg = MCpptree->GetPPSegment();
								double formerE = MCpptree->GetPPEnergy();
								std::cout << "Reading entry: "<<formerEvt <<", event No. " << MCpptree->GetPPEvent() << ", energy: " << formerE << " MeV, PSD: " << formerPSD << ", Segment: " << (formerSeg) << ", Y: " << std::abs(StartPos-formerPos) <<", Time: " << formerTime<< ", Time Diff: " << ((formerTime-StartTime)/1000) <<" \n";
								if (formerSeg == Segment && std::abs(StartPos-formerPos) < 50 && formerPSD > 0 && formerPSD < PSDGMax0 && formerE > ENoise && formerE<10 ) {
									hTime0->Fill((formerTime-StartTime)/1000);
									formerCount += 1;
								}
								
									formerEvt = formerEvt-1;
									MCpptree->GetEntry_PP_Tree(formerEvt);
									formerTime = MCpptree->GetPPTime();
							}
							std::cout<<"Former event number: "<< formerCount << std::endl;	

							MCpptree->GetEntry_PP_Tree(latterEvt);
							double latterTime = MCpptree->GetPPTime();
							//std::cout<<"Latter events:"<<std::endl;
							int latterCount = 0;
							while ((latterTime-StartTime)>0 && (latterTime-StartTime)<200000) {
								double latterPos = MCpptree->GetPPY();
								double latterPSD = MCpptree->GetPPPSD();
								int latterSeg = MCpptree->GetPPSegment();
								double latterE = MCpptree->GetPPEnergy();
								std::cout << "Reading entry: "<<latterEvt <<", event No. " << MCpptree->GetPPEvent() << ", energy: " << latterE << " MeV, PSD: " << latterPSD << ", Segment: " << (latterSeg) << ", Y: " << std::abs(StartPos-latterPos) <<", Time: " << latterTime<< ", Time Diff: " << ((latterTime-StartTime)/1000) <<" \n";
								if (latterSeg == Segment && std::abs(StartPos-latterPos) < 50 && latterPSD > 0 && latterPSD < PSDGMax0 && latterE > ENoise && latterE< 10)  {
									hTime0->Fill((latterTime-StartTime)/1000);
									latterCount += 1;
								}
									latterEvt = latterEvt+1;
									MCpptree->GetEntry_PP_Tree(latterEvt);
									latterTime = MCpptree->GetPPTime();
							}
							std::cout<<"Latter event number: "<< latterCount << std::endl;	
							
					//if (Segment == 1 && PSD > PSDGMax1) hELi1->Fill(Energy);
					}
				}
			}
	  }

		TFile* gInputFile = TFile::Open(InputFile.c_str());
		TVectorD* vRunTime = (TVectorD*)gInputFile->Get("runtime");
		PhysPulseRunTime += (*vRunTime)[0];
		delete MCpptree;
	}

	std::cout << "Total run time is "<< PhysPulseRunTime << "s.. \n";

	TFile* gOutputFile = TFile::Open(OutputFile.c_str(), "RECREATE");
	hEGammaTot->Sumw2();
	hEGamma0->Sumw2();
	hEGamma1->Sumw2();
	hELi0->Sumw2();
	hELi1->Sumw2();

	hEGammaTot->Scale(1./PhysPulseRunTime);
	hEGamma0->Scale(1./PhysPulseRunTime);
	hEGamma1->Scale(1./PhysPulseRunTime);
	hELi0->Scale(1./PhysPulseRunTime);
	hELi1->Scale(1./PhysPulseRunTime);

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
