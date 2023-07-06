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
  std::vector<float> fTrackVx;
  std::vector<float> fTrackVy;
  std::vector<float> fTrackVz;
  std::vector<float> fTrackPx;
  std::vector<float> fTrackPy;
  std::vector<float> fTrackPz;
  std::vector<int> fTrackHitsDetID;
  std::vector<int> fTrackHitsChan;
  std::vector<float> fTrackHitsPos;
  std::vector<float> fTrackHitsDrift;
  std::vector<int> fTrackHitsSign;
  std::vector<float> fTrackHitsTDC;
  std::vector<float> fTrackHitsResidual;
  std::vector<float> fTrackInvPTgt;

  Int_t fNDimuons;
  std::vector<float> fDimMass;
  std::vector<float> fDimPT;
  std::vector<float> fDimXF;
  std::vector<float> fDimX1;
  std::vector<float> fDimX2;
  std::vector<float> fDimCostheta;
  std::vector<float> fDimPhi;
  std::vector<float> fDimMassSingle;
  std::vector<float> fDimChi2Single;
  std::vector<float> fDimVx;
  std::vector<float> fDimVy;
  std::vector<float> fDimVz;
  std::vector<float> fDimPosVx;
  std::vector<float> fDimPosVy;
  std::vector<float> fDimPosVz;
  std::vector<float> fDimNegVx;
  std::vector<float> fDimNegVy;
  std::vector<float> fDimNegVz;
  std::vector<float> fDimPosE;
  std::vector<float> fDimPosPx;
  std::vector<float> fDimPosPy;
  std::vector<float> fDimPosPz;
  std::vector<float> fDimNegE;
  std::vector<float> fDimNegPx;
  std::vector<float> fDimNegPy;
  std::vector<float> fDimNegPz;
  std::vector<float> fDimPosSingleE;
  std::vector<float> fDimPosSinglePx;
  std::vector<float> fDimPosSinglePy;
  std::vector<float> fDimPosSinglePz;
  std::vector<float> fDimNegSingleE;
  std::vector<float> fDimNegSinglePx;
  std::vector<float> fDimNegSinglePy;
  std::vector<float> fDimNegSinglePz;
  std::vector<float> fDimChi2Vtx;
  std::vector<float> fDimChi2KF;
  
  ClassDef(ORoutput_tree, 1)
}; 

