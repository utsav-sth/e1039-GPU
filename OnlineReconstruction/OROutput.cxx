#include <iostream>
#include <vector>
#include <TChain.h>
#include <TROOT.h>
#include "OROutput.h"

ClassImp(ORoutput_tree)

ORoutput_tree::ORoutput_tree(const char* filename)
{
  fFile = new TFile(filename, "RECREATE");
  fTree =  new TTree("tree", "Online reconstruction output");
  Init();
}

ORoutput_tree::~ORoutput_tree()
{
  
}

void ORoutput_tree::Init()
{
  Clear();
  fTree->Branch("evtid", &(fEventID));
  
  fTree->Branch("nhits", &(fNHits));
  fTree->Branch("hit.detid", &(fHitDetID));
  fTree->Branch("hit.chan", &(fHitChan));
  fTree->Branch("hit.pos", &(fHitPos));
  fTree->Branch("hit.tdc", &(fHitTDC));
  fTree->Branch("hit.drift", &(fHitDrift));
  
  fTree->Branch("ntracks", &(fNTracks));
  fTree->Branch("track.stid", &(fTrackStID));
  fTree->Branch("track.nhits", &(fTrackNHits));
  fTree->Branch("track.chi2", &(fTrackChi2));
  fTree->Branch("track.x0", &(fTrackX0));
  fTree->Branch("track.y0", &(fTrackY0));
  fTree->Branch("track.tx", &(fTrackTX));
  fTree->Branch("track.ty", &(fTrackTY));
  fTree->Branch("track.invp", &(fTrackInvP));
  fTree->Branch("track.err_x0", &(fTrackErrX0));
  fTree->Branch("track.err_y0", &(fTrackErrY0));
  fTree->Branch("track.err_tx", &(fTrackErrTX));
  fTree->Branch("track.err_ty", &(fTrackErrTY));
  fTree->Branch("track.err_invp", &(fTrackErrInvP));
  fTree->Branch("track.charge", &(fTrackCharge));
  fTree->Branch("track.hitdetid", &(fTrackHitsDetID));
  fTree->Branch("track.hitchan", &(fTrackHitsChan));
  fTree->Branch("track.hitpos", &(fTrackHitsPos));
  fTree->Branch("track.hitdrift", &(fTrackHitsDrift));
  fTree->Branch("track.hitsign", &(fTrackHitsSign));
  fTree->Branch("track.vx", &(fTrackVx));
  fTree->Branch("track.vy", &(fTrackVy));
  fTree->Branch("track.vz", &(fTrackVz));
  fTree->Branch("track.px", &(fTrackPx));
  fTree->Branch("track.py", &(fTrackPy));
  fTree->Branch("track.pz", &(fTrackPz));
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
  fTrackHitsDrift.clear();
  fTrackHitsSign.clear();
  fTrackVx.clear();
  fTrackVy.clear();
  fTrackVz.clear();
  fTrackPx.clear();
  fTrackPy.clear();
  fTrackPz.clear();
}

void ORoutput_tree::FillTree()
{
  fTree->Fill();
}

void ORoutput_tree::Close()
{
  fTree->Write();
  fFile->Close();
}

