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
  fTree->Branch("track.vx", &(fTrackVx));
  fTree->Branch("track.vy", &(fTrackVy));
  fTree->Branch("track.vz", &(fTrackVz));
  fTree->Branch("track.px", &(fTrackPx));
  fTree->Branch("track.py", &(fTrackPy));
  fTree->Branch("track.pz", &(fTrackPz));
  fTree->Branch("track.hitdetid", &(fTrackHitsDetID));
  fTree->Branch("track.hitchan", &(fTrackHitsChan));
  fTree->Branch("track.hitpos", &(fTrackHitsPos));
  fTree->Branch("track.hitdrift", &(fTrackHitsDrift));
  fTree->Branch("track.hitsign", &(fTrackHitsSign));
  fTree->Branch("track.hittdc", &(fTrackHitsTDC));
  fTree->Branch("track.hitresid", &(fTrackHitsResidual));
  fTree->Branch("track.invp_tgt", &(fTrackInvPTgt));
  
  fTree->Branch("ndimuons", &(fNDimuons));
  fTree->Branch("dimuon.mass", &(fDimMass));
  fTree->Branch("dimuon.pT", &(fDimPT)); 
  fTree->Branch("dimuon.xF", &(fDimXF));
  fTree->Branch("dimuon.x1", &(fDimX1));
  fTree->Branch("dimuon.x2", &(fDimX2));
  fTree->Branch("dimuon.costheta", &(fDimCostheta));
  fTree->Branch("dimuon.phi", &(fDimPhi));
  fTree->Branch("dimuon.mass_single", &(fDimMassSingle));
  fTree->Branch("dimuon.chi2_single", &(fDimChi2Single));
  fTree->Branch("dimuon.vx", &(fDimVx));
  fTree->Branch("dimuon.vy", &(fDimVy));
  fTree->Branch("dimuon.vz", &(fDimVz));
  fTree->Branch("dimuon.pos_vx", &(fDimPosVx));
  fTree->Branch("dimuon.pos_vy", &(fDimPosVy));
  fTree->Branch("dimuon.pos_vz", &(fDimPosVz));
  fTree->Branch("dimuon.neg_vx", &(fDimNegVx));
  fTree->Branch("dimuon.neg_vy", &(fDimNegVy));
  fTree->Branch("dimuon.neg_vz", &(fDimNegVz));
  fTree->Branch("dimuon.pos_e", &(fDimPosE));
  fTree->Branch("dimuon.pos_px", &(fDimPosPx));
  fTree->Branch("dimuon.pos_py", &(fDimPosPy));
  fTree->Branch("dimuon.pos_pz", &(fDimPosPz));
  fTree->Branch("dimuon.neg_e", &(fDimNegE));
  fTree->Branch("dimuon.neg_px", &(fDimNegPx));
  fTree->Branch("dimuon.neg_py", &(fDimNegPy));
  fTree->Branch("dimuon.neg_pz", &(fDimNegPz));
  fTree->Branch("dimuon.pos_single_e", &(fDimPosSingleE));
  fTree->Branch("dimuon.pos_single_px", &(fDimPosSinglePx));
  fTree->Branch("dimuon.pos_single_py", &(fDimPosSinglePy));
  fTree->Branch("dimuon.pos_single_pz", &(fDimPosSinglePz));
  fTree->Branch("dimuon.neg_single_e", &(fDimNegSingleE));
  fTree->Branch("dimuon.neg_single_px", &(fDimNegSinglePx));
  fTree->Branch("dimuon.neg_single_py", &(fDimNegSinglePy));
  fTree->Branch("dimuon.neg_single_pz", &(fDimNegSinglePz));
  fTree->Branch("dimuon.chi2_vtx", &(fDimChi2Vtx));
  fTree->Branch("dimuon.chi2_kf", &(fDimChi2KF));
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
  fTrackVx.clear();
  fTrackVy.clear();
  fTrackVz.clear();
  fTrackPx.clear();
  fTrackPy.clear();
  fTrackPz.clear();
  fTrackHitsDetID.clear();
  fTrackHitsChan.clear();
  fTrackHitsPos.clear();
  fTrackHitsDrift.clear();
  fTrackHitsSign.clear();
  fTrackHitsTDC.clear();
  fTrackHitsResidual.clear();
  fTrackInvPTgt.clear();
  
  fNDimuons = 0;
  fDimMass.clear();
  fDimPT.clear();
  fDimXF.clear();
  fDimX1.clear();
  fDimX2.clear();
  fDimCostheta.clear();
  fDimPhi.clear();
  fDimMassSingle.clear();
  fDimChi2Single.clear();
  fDimVx.clear();
  fDimVy.clear();
  fDimVz.clear();
  fDimPosVx.clear();
  fDimPosVy.clear();
  fDimPosVz.clear();
  fDimNegVx.clear();
  fDimNegVy.clear();
  fDimNegVz.clear();
  fDimPosE.clear();
  fDimPosPx.clear();
  fDimPosPy.clear();
  fDimPosPz.clear();
  fDimNegE.clear();
  fDimNegPx.clear();
  fDimNegPy.clear();
  fDimNegPz.clear();
  fDimPosSingleE.clear();
  fDimPosSinglePx.clear();
  fDimPosSinglePy.clear();
  fDimPosSinglePz.clear();
  fDimNegSingleE.clear();
  fDimNegSinglePx.clear();
  fDimNegSinglePy.clear();
  fDimNegSinglePz.clear();
  fDimChi2Vtx.clear();
  fDimChi2KF.clear();
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

