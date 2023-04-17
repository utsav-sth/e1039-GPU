#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include "LoadInput.h"

class ORoutput_tree
{
public:
  //ORoutput_tree();
  ORoutput_tree(const char* filename);
  ~ORoutput_tree();

  void Init();
  void Clear();
  void FillTree();
  void Write();
  void Close();
  
  TFile* fFile;
  TTree* fTree;

  Int_t fEventID;
  Int_t fNHits;
  std::vector<int> fHitDetID;
  std::vector<int> fHitChan;
  std::vector<float> fHitPos;
  std::vector<float> fHitTDC;
  std::vector<float> fHitDrift;

  Int_t fNTracks;
  std::vector<int> fTrackStID;
  std::vector<int> fTrackNHits;
  std::vector<float> fTrackChi2;
  std::vector<float> fTrackX0;
  std::vector<float> fTrackY0;
  std::vector<float> fTrackTX;
  std::vector<float> fTrackTY;
  std::vector<float> fTrackInvP;
  std::vector<float> fTrackErrX0;
  std::vector<float> fTrackErrY0;
  std::vector<float> fTrackErrTX;
  std::vector<float> fTrackErrTY;
  std::vector<float> fTrackErrInvP;
  std::vector<int> fTrackCharge;
  std::vector<int> fTrackHitsDetID;
  std::vector<int> fTrackHitsChan;
  std::vector<float> fTrackHitsPos;
  std::vector<float> fTrackHitsDrift;
  std::vector<float> fTrackVx;
  std::vector<float> fTrackVy;
  std::vector<float> fTrackVz;
  std::vector<float> fTrackPx;
  std::vector<float> fTrackPy;
  std::vector<float> fTrackPz;

  ClassDef(ORoutput_tree, 1)
}; 

