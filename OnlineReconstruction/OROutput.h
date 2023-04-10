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
  std::vector<float> fTrackHitsTDC;
  std::vector<float> fTrackHitsDrift;
  std::vector<float> fTrackVx;
  std::vector<float> fTrackVy;
  std::vector<float> fTrackVz;
  std::vector<float> fTrackPx;
  std::vector<float> fTrackPy;
  std::vector<float> fTrackPz;

  /*
  TBranch* b_NHits;
  TBranch* b_HitDetID;
  TBranch* b_HitChan;
  TBranch* b_HitPos;
  TBranch* b_HitTDC;
  TBranch* b_HitDrift;
  
  TBranch* b_NTracks;
  TBranch* b_TrackStID;
  TBranch* b_TrackNHits;
  TBranch* b_TrackChi2;
  TBranch* b_TrackX0;
  TBranch* b_TrackY0;
  TBranch* b_TrackTX;
  TBranch* b_TrackTY;
  TBranch* b_TrackInvP;
  TBranch* b_TrackErrX0;
  TBranch* b_TrackErrY0;
  TBranch* b_TrackErrTX;
  TBranch* b_TrackErrTY;
  TBranch* b_TrackErrInvP;
  TBranch* b_TrackCharge;
  TBranch* b_TrackHitsDetID;
  TBranch* b_TrackHitsChan;
  TBranch* b_TrackHitsPos;
  TBranch* b_TrackHitsTDC;
  TBranch* b_TrackHitsDrift;
  TBranch* b_TrackVx;
  TBranch* b_TrackVy;
  TBranch* b_TrackVz;
  TBranch* b_TrackPx;
  TBranch* b_TrackPy;
  TBranch* b_TrackPz;
  */
  
  ClassDef(ORoutput_tree, 1)
}; 

