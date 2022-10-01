#include <iostream>
#include <vector>
#include <TChain.h>
#include <TBranch.h>
#include <TROOT.h>
#include "LoadInput.h"

class ORoutput_tree
{
public:
  ORoutput_tree();
  ~ORoutput_tree();

  void Init();
  void Clear();
  void Write();
  
  TChain* fChain;
  Int_t fNhitsReduced[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  //std::vector<Hit> fHitsReduced;
  
// private:
//   TBranch* b_nhitsreduced;
//   TBranch* b_hitsreduced;
  
  ClassDef(ORoutput_tree, 1)
}; 

