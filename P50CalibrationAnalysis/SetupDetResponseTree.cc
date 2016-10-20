#include "SetupDetResponseTree.hh"

PhysPulseTree::PhysPulseTree(){	
	gPhysPulse_Tree = new TChain("PhysPulse");
}

void PhysPulseTree::SetupPhysPulseTree(std::string fInputFile){
	
	gPhysPulse_Tree = new TChain("PhysPulse");
	gPhysPulse_Tree->AddFile(fInputFile.c_str());

	gPhysPulse_Tree->SetBranchAddress("evt", &gPhysPulse_Event);
	gPhysPulse_Tree->SetBranchAddress("seg", &gPhysPulse_Segment);
	gPhysPulse_Tree->SetBranchAddress("E", &gPhysPulse_Energy);
	gPhysPulse_Tree->SetBranchAddress("t", &gPhysPulse_Time);
	gPhysPulse_Tree->SetBranchAddress("dt", &gPhysPulse_dT);
	gPhysPulse_Tree->SetBranchAddress("PE", &gPhysPulse_PE);
	gPhysPulse_Tree->SetBranchAddress("y", &gPhysPulse_Y);
	gPhysPulse_Tree->SetBranchAddress("PSD", &gPhysPulse_PSD);
	gPhysPulse_Tree->SetBranchAddress("PID", &gPhysPulse_PID);

	gPhysPulse_Tree->GetEntry(0);
}

int PhysPulseTree::GetEntries_PP_Tree(){
	return gPhysPulse_Tree->GetEntries();
}

int PhysPulseTree::GetEntry_PP_Tree(int ev){
	return gPhysPulse_Tree->GetEntry(ev);
}

ULong64_t PhysPulseTree::GetPPEvent(){
	return gPhysPulse_Event;
}

Int_t PhysPulseTree::GetPPSegment(){
	return gPhysPulse_Segment;
}

Float_t PhysPulseTree::GetPPEnergy(){
	return gPhysPulse_Energy;
}

Double_t PhysPulseTree::GetPPTime(){
	return gPhysPulse_Time;
}

Double_t PhysPulseTree::GetPPdT(){
	return gPhysPulse_dT;
}

Double_t PhysPulseTree::GetPPPE(int i){
	return gPhysPulse_PE[i];
}

Double_t PhysPulseTree::GetPPY(){
	return gPhysPulse_Y;
}

Double_t PhysPulseTree::GetPPPSD(){
	return gPhysPulse_PSD;
}

Int_t PhysPulseTree::GetPPPID(){
	return gPhysPulse_PID;
}
