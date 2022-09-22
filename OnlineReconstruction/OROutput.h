#include <iostream>
#include <vector>
#include <TTree.h>
#include <TROOT.h>
#include "LoadInput.h"

class ORoutput_tree: public TTree
{
public:
    ORoutput_tree();
    ~ORoutput_tree();

  void Init();
  
    Int_t fNhitsReduced[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
    std::vector<Hit> fHitsReduced;
    
    ClassDef(ORoutput_tree, 1)
}; 

