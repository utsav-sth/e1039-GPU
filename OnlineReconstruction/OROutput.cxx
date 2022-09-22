#include <iostream>
#include <vector>
#include <TTree.h>
#include <TROOT.h>
#include "Output.h"

ClassImp(ORoutput_tree)

ORoutput_tree::ORoutput_tree() : 
{
   for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
     {
       fNHitsReduced[i] = 0;
     }
}

ORoutput_tree::ORoutput_tree~()
{
  
}

void ORoutput_tree::Init()
{
  SetBranchAddress(fNhitsReduced,"Nreducedhits[55]");
  SetBranchAddress(fHitsReduced,"reducedhits");
}
