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
  
  fChain->SetBranchAddress("nhits", &fNHits);
  fChain->SetBranchAddress("hit.detid", &fHitDetID);
  fChain->SetBranchAddress("hit.chan", &fHitChan);
  fChain->SetBranchAddress("hit.pos", &fHitPos);
  fChain->SetBranchAddress("hit.tdc", &fHitTDC);
  fChain->SetBranchAddress("hit.drift", &fHitDrift);
  
  fChain->SetBranchAddress("ntracks", &fNTracks);
  fChain->SetBranchAddress("track.stid", &fTrackStID);
  fChain->SetBranchAddress("track.nhits", &fTrackNHits);
  fChain->SetBranchAddress("track.chi2", &fTrackChi2);
  fChain->SetBranchAddress("track.x0", &fTrackX0);
  fChain->SetBranchAddress("track.y0", &fTrackY0);
  fChain->SetBranchAddress("track.tx", &fTrackTX);
  fChain->SetBranchAddress("track.ty", &fTrackTY);
  fChain->SetBranchAddress("track.invp", &fTrackInvP);
  fChain->SetBranchAddress("track.err_x0", &fTrackErrX0);
  fChain->SetBranchAddress("track.err_y0", &fTrackErrY0);
  fChain->SetBranchAddress("track.err_tx", &fTrackErrTX);
  fChain->SetBranchAddress("track.err_ty", &fTrackErrTY);
  fChain->SetBranchAddress("track.err_invp", &fTrackErrInvP);
  fChain->SetBranchAddress("track.charge", &fTrackCharge);
  fChain->SetBranchAddress("track.hitdetid", &fTrackHitsDetID);
  fChain->SetBranchAddress("track.hitchan", &fTrackHitsChan);
  fChain->SetBranchAddress("track.hitpos", &fTrackHitsPos);
  fChain->SetBranchAddress("track.hittdc", &fTrackHitsTDC);
  fChain->SetBranchAddress("track.hitdrift", &fTrackHitsDrift);
  fChain->SetBranchAddress("track.vx", &fTrackVx);
  fChain->SetBranchAddress("track.vy", &fTrackVy);
  fChain->SetBranchAddress("track.vz", &fTrackVz);
  fChain->SetBranchAddress("track.px", &fTrackPx);
  fChain->SetBranchAddress("track.py", &fTrackPy);
  fChain->SetBranchAddress("track.pz", &fTrackPz);
}

void ORoutput_tree::Clear()
{
  fNHits = 0;
  fHitDetID.clear();
  fHitChan.clear();
  fHitPos.clear();
  fHitTDC.clear();
  fHitDrift.clear();
  fNTracks = 0;
  fTrackStID.clear();
  fTrackNHits.clear();
  fTrackChi2.clear();
  fTrackX0.clear();
  fTrackY0.clear();
  fTrackTX.clear();
  fTrackTY.clear();
  fTrackInvP.clear();
  fTrackErrX0.clear();
  fTrackErrY0.clear();
  fTrackErrTX.clear();
  fTrackErrTY.clear();
  fTrackErrInvP.clear();
  fTrackCharge.clear();
  fTrackHitsDetID.clear();
  fTrackHitsChan.clear();
  fTrackHitsPos.clear();
  fTrackHitsTDC.clear();
  fTrackHitsDrift.clear();
}

void ORoutput_tree::Write()
{
  fChain->Write("");
}
