#ifndef DETRESPONSE_HH
#define DETRESPONSE_HH

#include <iostream>
#include "TH1F.h"

#include "Event.hh"

class DetResponse{

public:
  DetResponse();
  void FillPrimEvent(ParticleVertex*, Long64_t);
  void ClearEnergy();
  void FillHisto();
  void FillDetEvent(IoniCluster*, Long64_t);
  void DrawHisto();
	double energy0Out = 0;
	double energy1Out = 0;

public:
  TH1F* E_prim;
  TH1F* E_dep_tot;
	TH1F* E_dep_0;
	TH1F* E_dep_1;

protected:
  int testEvt = 0;
	int testEvtIoni = 0;
	double primE = 0;
  double energy[2] = {0,0};
};

#endif
