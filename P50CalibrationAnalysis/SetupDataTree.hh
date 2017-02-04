#ifndef SETUPDATATREE_HH
#define SETUPDATATREE_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TVectorD.h"

class PhysPulseTree {
public:
	PhysPulseTree(std::string);

	int GetEntries_PP_Tree();

	int GetEntry_PP_Tree(int);

	Long64_t GetPPEvent();

	Int_t GetPPSegment();

	Float_t GetPPEnergy();

	Double_t GetPPTime();

	Float_t GetPPdT();

	Double_t GetPPPE(int);

	Float_t GetPPY();

	Float_t GetPPPSD();

	Int_t GetPPPID();
	
	Double_t gPhysPulse_RunTime;
	
	Double_t gPhysPulse_AbsTime;

protected: 
	TChain* gPhysPulse_Tree = NULL;

	Long64_t gPhysPulse_Event = 0;

	Int_t gPhysPulse_Segment = 0;

	Float_t gPhysPulse_Energy = 0;

	Double_t gPhysPulse_Time = 0;

	Float_t gPhysPulse_dT = 0;

	Double_t gPhysPulse_PE[2] = {0.};

	Float_t gPhysPulse_Y = 0;

	Float_t gPhysPulse_PSD = 0;

	Int_t gPhysPulse_PID = 0;
};

#endif
