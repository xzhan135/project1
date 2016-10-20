#ifndef SETUPDETRESPONSETREE_HH
#define SETUPDETRESPONSETREE_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TString.h"

class PhysPulseTree {

public:
	PhysPulseTree();

	void SetupPhysPulseTree(std::string);

	int GetEntries_PP_Tree();

	int GetEntry_PP_Tree(int);

	ULong64_t GetPPEvent();

	Int_t GetPPSegment();

	Float_t GetPPEnergy();

	Double_t GetPPTime();

	Double_t GetPPdT();

	Double_t GetPPPE(int);

	Double_t GetPPY();

	Double_t GetPPPSD();

	Int_t GetPPPID();

protected: 
	TChain* gPhysPulse_Tree = NULL;

	ULong64_t gPhysPulse_Event = 0;

	Int_t gPhysPulse_Segment = 0;

	Double_t gPhysPulse_Energy = 0;

	Double_t gPhysPulse_Time = 0;

	Double_t gPhysPulse_dT = 0;

	Double_t gPhysPulse_PE[2] = {0.};

	Double_t gPhysPulse_Y = 0;

	Double_t gPhysPulse_PSD = 0;

	Int_t gPhysPulse_PID = 0;
};

#endif
