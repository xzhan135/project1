#include "SetupDataTree.hh"

PhysPulseTree::PhysPulseTree(std::string inputFile){
	gPhysPulse_Tree=new TChain("PhysPulse");

	std::ifstream inputList(fileList.c_str());

	gPhysPulse_Tree->Add(inputFile.c_str());

	gPhysPulse_Tree->SetBranchAddress("evt", &gPhysPulse_Event);
	gPhysPulse_Tree->SetBranchAddress("seg", &gPhysPulse_Segment);
	gPhysPulse_Tree->SetBranchAddress("E", &gPhysPulse_Energy);
	gPhysPulse_Tree->SetBranchAddress("t", &gPhysPulse_Time);
	gPhysPulse_Tree->SetBranchAddress("dt", &gPhysPulse_dT);
	gPhysPulse_Tree->SetBranchAddress("PE", &gPhysPulse_PE);
	gPhysPulse_Tree->SetBranchAddress("y", &gPhysPulse_Y);
	gPhysPulse_Tree->SetBranchAddress("PSD", &gPhysPulse_PSD);
	gPhysPulse_Tree->SetBranchAddress("PID", &gPhysPulse_PID);
}

int PhysPulseTree::GetEntries_PP_Tree(){
	return gPhysPulse_Tree->GetEntries();
}

int PhysPulseTree::GetEntry_PP_Tree(int ev){
	return gPhysPulse_Tree->GetEntry(ev);
}

Long64_t PhysPulseTree::GetPPEvent(){
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

Float_t PhysPulseTree::GetPPdT(){
	return gPhysPulse_dT;
}

Double_t PhysPulseTree::GetPPPE(int i){
	return gPhysPulse_PE[i];
}

Float_t PhysPulseTree::GetPPY(){
	return gPhysPulse_Y;
}

Float_t PhysPulseTree::GetPPPSD(){
	return gPhysPulse_PSD;
}

Int_t PhysPulseTree::GetPPPID(){
	return gPhysPulse_PID;
}
