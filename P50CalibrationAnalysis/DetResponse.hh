#ifndef DETRESPONSE_HH
#define DETRESPONSE_HH

#include "TH1F.h"

#include "Event.hh"

class DetResponse{

public:
  DetResponse();
  void FillPrimEvent(ParticleVertex*, Long64_t);
  void FillDetEvent(IoniCluster*, Long64_t);
  void DrawHisto();

public:
  TH1F* E_prim;
  TH1F* E_dep_tot;
	TH1F* E_dep_0;
	TH1F* E_dep_1;

protected:
  int testEvt = 0;
  double energy = 0;
};

#endif
