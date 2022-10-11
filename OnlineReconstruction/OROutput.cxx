#include <iostream>
#include <vector>
#include <TChain.h>
#include <TROOT.h>
#include "OROutput.h"

ClassImp(ORoutput_tree)

ORoutput_tree::ORoutput_tree()
{
  fChain =  new TChain("ORoutput", "ORout");
  Init();
}

ORoutput_tree::~ORoutput_tree()
{
  
}

void ORoutput_tree::Init()
{
  Clear();
  //fHitsReduced.clear();
  
  fChain->SetBranchAddress("Nreducedhits[55]", fNhitsReduced);
  //fChain->SetBranchAddress("reducedhits.", fHitsReduced.);
			   
}

void ORoutput_tree::Clear()
{
  for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    {
      fNhitsReduced[i] = 0;
    }
}

void ORoutput_tree::Write()
{
  fChain->Write("");
  //fChain->SetBranchAddress("reducedhits.", fHitsReduced.);
			   
}
