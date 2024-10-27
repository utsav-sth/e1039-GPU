//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 10 23:14:43 2022 by ROOT version 6.16/00
// from TTree T/titled by PHOOL
// found on file: DST.root
//////////////////////////////////////////////////////////

#ifndef e1039data_tree_h
#define e1039data_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class e1039data_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxDST_SRecEvent_fAllTracks = 2;
   static constexpr Int_t kMaxDST_SRecEvent_fDimuons = 1;

   // Declaration of leaf types
 //PHHepMCGenEventMap *DST_PHHepMCGenEventMap;
   UInt_t          DST_PHHepMCGenEventMap_fUniqueID;
   UInt_t          DST_PHHepMCGenEventMap_fBits;
   map<int,PHHepMCGenEvent*> DST_PHHepMCGenEventMap__map;
 //PHG4InEvent     *DST_PHG4INEVENT;
   UInt_t          DST_PHG4INEVENT_fUniqueID;
   UInt_t          DST_PHG4INEVENT_fBits;
   map<int,PHG4VtxPoint*> DST_PHG4INEVENT_vtxlist;
   multimap<int,PHG4Particle*> DST_PHG4INEVENT_particlelist;
   map<PHG4Particle*,int> DST_PHG4INEVENT_embedded_particlelist;
 //PHG4HitContainer *DST_G4HIT_D0U;
   UInt_t          DST_G4HIT_D0U_fUniqueID;
   UInt_t          DST_G4HIT_D0U_fBits;
   Int_t           DST_G4HIT_D0U_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0U_hitmap;
   set<unsigned int> DST_G4HIT_D0U_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0U_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D0Up;
   UInt_t          DST_G4HIT_D0Up_fUniqueID;
   UInt_t          DST_G4HIT_D0Up_fBits;
   Int_t           DST_G4HIT_D0Up_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0Up_hitmap;
   set<unsigned int> DST_G4HIT_D0Up_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0Up_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D0X;
   UInt_t          DST_G4HIT_D0X_fUniqueID;
   UInt_t          DST_G4HIT_D0X_fBits;
   Int_t           DST_G4HIT_D0X_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0X_hitmap;
   set<unsigned int> DST_G4HIT_D0X_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0X_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D0Xp;
   UInt_t          DST_G4HIT_D0Xp_fUniqueID;
   UInt_t          DST_G4HIT_D0Xp_fBits;
   Int_t           DST_G4HIT_D0Xp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0Xp_hitmap;
   set<unsigned int> DST_G4HIT_D0Xp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0Xp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D0V;
   UInt_t          DST_G4HIT_D0V_fUniqueID;
   UInt_t          DST_G4HIT_D0V_fBits;
   Int_t           DST_G4HIT_D0V_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0V_hitmap;
   set<unsigned int> DST_G4HIT_D0V_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0V_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D0Vp;
   UInt_t          DST_G4HIT_D0Vp_fUniqueID;
   UInt_t          DST_G4HIT_D0Vp_fBits;
   Int_t           DST_G4HIT_D0Vp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D0Vp_hitmap;
   set<unsigned int> DST_G4HIT_D0Vp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D0Vp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2V;
   UInt_t          DST_G4HIT_D2V_fUniqueID;
   UInt_t          DST_G4HIT_D2V_fBits;
   Int_t           DST_G4HIT_D2V_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2V_hitmap;
   set<unsigned int> DST_G4HIT_D2V_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2V_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2Vp;
   UInt_t          DST_G4HIT_D2Vp_fUniqueID;
   UInt_t          DST_G4HIT_D2Vp_fBits;
   Int_t           DST_G4HIT_D2Vp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2Vp_hitmap;
   set<unsigned int> DST_G4HIT_D2Vp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2Vp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2Xp;
   UInt_t          DST_G4HIT_D2Xp_fUniqueID;
   UInt_t          DST_G4HIT_D2Xp_fBits;
   Int_t           DST_G4HIT_D2Xp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2Xp_hitmap;
   set<unsigned int> DST_G4HIT_D2Xp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2Xp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2X;
   UInt_t          DST_G4HIT_D2X_fUniqueID;
   UInt_t          DST_G4HIT_D2X_fBits;
   Int_t           DST_G4HIT_D2X_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2X_hitmap;
   set<unsigned int> DST_G4HIT_D2X_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2X_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2U;
   UInt_t          DST_G4HIT_D2U_fUniqueID;
   UInt_t          DST_G4HIT_D2U_fBits;
   Int_t           DST_G4HIT_D2U_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2U_hitmap;
   set<unsigned int> DST_G4HIT_D2U_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2U_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D2Up;
   UInt_t          DST_G4HIT_D2Up_fUniqueID;
   UInt_t          DST_G4HIT_D2Up_fBits;
   Int_t           DST_G4HIT_D2Up_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D2Up_hitmap;
   set<unsigned int> DST_G4HIT_D2Up_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D2Up_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pVp;
   UInt_t          DST_G4HIT_D3pVp_fUniqueID;
   UInt_t          DST_G4HIT_D3pVp_fBits;
   Int_t           DST_G4HIT_D3pVp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pVp_hitmap;
   set<unsigned int> DST_G4HIT_D3pVp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pVp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pV;
   UInt_t          DST_G4HIT_D3pV_fUniqueID;
   UInt_t          DST_G4HIT_D3pV_fBits;
   Int_t           DST_G4HIT_D3pV_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pV_hitmap;
   set<unsigned int> DST_G4HIT_D3pV_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pV_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pXp;
   UInt_t          DST_G4HIT_D3pXp_fUniqueID;
   UInt_t          DST_G4HIT_D3pXp_fBits;
   Int_t           DST_G4HIT_D3pXp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pXp_hitmap;
   set<unsigned int> DST_G4HIT_D3pXp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pXp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pX;
   UInt_t          DST_G4HIT_D3pX_fUniqueID;
   UInt_t          DST_G4HIT_D3pX_fBits;
   Int_t           DST_G4HIT_D3pX_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pX_hitmap;
   set<unsigned int> DST_G4HIT_D3pX_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pX_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pUp;
   UInt_t          DST_G4HIT_D3pUp_fUniqueID;
   UInt_t          DST_G4HIT_D3pUp_fBits;
   Int_t           DST_G4HIT_D3pUp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pUp_hitmap;
   set<unsigned int> DST_G4HIT_D3pUp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pUp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3pU;
   UInt_t          DST_G4HIT_D3pU_fUniqueID;
   UInt_t          DST_G4HIT_D3pU_fBits;
   Int_t           DST_G4HIT_D3pU_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3pU_hitmap;
   set<unsigned int> DST_G4HIT_D3pU_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3pU_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mVp;
   UInt_t          DST_G4HIT_D3mVp_fUniqueID;
   UInt_t          DST_G4HIT_D3mVp_fBits;
   Int_t           DST_G4HIT_D3mVp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mVp_hitmap;
   set<unsigned int> DST_G4HIT_D3mVp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mVp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mV;
   UInt_t          DST_G4HIT_D3mV_fUniqueID;
   UInt_t          DST_G4HIT_D3mV_fBits;
   Int_t           DST_G4HIT_D3mV_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mV_hitmap;
   set<unsigned int> DST_G4HIT_D3mV_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mV_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mXp;
   UInt_t          DST_G4HIT_D3mXp_fUniqueID;
   UInt_t          DST_G4HIT_D3mXp_fBits;
   Int_t           DST_G4HIT_D3mXp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mXp_hitmap;
   set<unsigned int> DST_G4HIT_D3mXp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mXp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mX;
   UInt_t          DST_G4HIT_D3mX_fUniqueID;
   UInt_t          DST_G4HIT_D3mX_fBits;
   Int_t           DST_G4HIT_D3mX_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mX_hitmap;
   set<unsigned int> DST_G4HIT_D3mX_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mX_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mUp;
   UInt_t          DST_G4HIT_D3mUp_fUniqueID;
   UInt_t          DST_G4HIT_D3mUp_fBits;
   Int_t           DST_G4HIT_D3mUp_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mUp_hitmap;
   set<unsigned int> DST_G4HIT_D3mUp_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mUp_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_D3mU;
   UInt_t          DST_G4HIT_D3mU_fUniqueID;
   UInt_t          DST_G4HIT_D3mU_fBits;
   Int_t           DST_G4HIT_D3mU_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_D3mU_hitmap;
   set<unsigned int> DST_G4HIT_D3mU_layers;
   map<unsigned int,unsigned int> DST_G4HIT_D3mU_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H1B;
   UInt_t          DST_G4HIT_H1B_fUniqueID;
   UInt_t          DST_G4HIT_H1B_fBits;
   Int_t           DST_G4HIT_H1B_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H1B_hitmap;
   set<unsigned int> DST_G4HIT_H1B_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H1B_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H1T;
   UInt_t          DST_G4HIT_H1T_fUniqueID;
   UInt_t          DST_G4HIT_H1T_fBits;
   Int_t           DST_G4HIT_H1T_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H1T_hitmap;
   set<unsigned int> DST_G4HIT_H1T_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H1T_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H1L;
   UInt_t          DST_G4HIT_H1L_fUniqueID;
   UInt_t          DST_G4HIT_H1L_fBits;
   Int_t           DST_G4HIT_H1L_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H1L_hitmap;
   set<unsigned int> DST_G4HIT_H1L_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H1L_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H1R;
   UInt_t          DST_G4HIT_H1R_fUniqueID;
   UInt_t          DST_G4HIT_H1R_fBits;
   Int_t           DST_G4HIT_H1R_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H1R_hitmap;
   set<unsigned int> DST_G4HIT_H1R_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H1R_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H2L;
   UInt_t          DST_G4HIT_H2L_fUniqueID;
   UInt_t          DST_G4HIT_H2L_fBits;
   Int_t           DST_G4HIT_H2L_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H2L_hitmap;
   set<unsigned int> DST_G4HIT_H2L_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H2L_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H2R;
   UInt_t          DST_G4HIT_H2R_fUniqueID;
   UInt_t          DST_G4HIT_H2R_fBits;
   Int_t           DST_G4HIT_H2R_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H2R_hitmap;
   set<unsigned int> DST_G4HIT_H2R_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H2R_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H2B;
   UInt_t          DST_G4HIT_H2B_fUniqueID;
   UInt_t          DST_G4HIT_H2B_fBits;
   Int_t           DST_G4HIT_H2B_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H2B_hitmap;
   set<unsigned int> DST_G4HIT_H2B_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H2B_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H2T;
   UInt_t          DST_G4HIT_H2T_fUniqueID;
   UInt_t          DST_G4HIT_H2T_fBits;
   Int_t           DST_G4HIT_H2T_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H2T_hitmap;
   set<unsigned int> DST_G4HIT_H2T_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H2T_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H3B;
   UInt_t          DST_G4HIT_H3B_fUniqueID;
   UInt_t          DST_G4HIT_H3B_fBits;
   Int_t           DST_G4HIT_H3B_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H3B_hitmap;
   set<unsigned int> DST_G4HIT_H3B_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H3B_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H3T;
   UInt_t          DST_G4HIT_H3T_fUniqueID;
   UInt_t          DST_G4HIT_H3T_fBits;
   Int_t           DST_G4HIT_H3T_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H3T_hitmap;
   set<unsigned int> DST_G4HIT_H3T_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H3T_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4Y1L;
   UInt_t          DST_G4HIT_H4Y1L_fUniqueID;
   UInt_t          DST_G4HIT_H4Y1L_fBits;
   Int_t           DST_G4HIT_H4Y1L_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4Y1L_hitmap;
   set<unsigned int> DST_G4HIT_H4Y1L_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4Y1L_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4Y1R;
   UInt_t          DST_G4HIT_H4Y1R_fUniqueID;
   UInt_t          DST_G4HIT_H4Y1R_fBits;
   Int_t           DST_G4HIT_H4Y1R_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4Y1R_hitmap;
   set<unsigned int> DST_G4HIT_H4Y1R_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4Y1R_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4Y2L;
   UInt_t          DST_G4HIT_H4Y2L_fUniqueID;
   UInt_t          DST_G4HIT_H4Y2L_fBits;
   Int_t           DST_G4HIT_H4Y2L_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4Y2L_hitmap;
   set<unsigned int> DST_G4HIT_H4Y2L_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4Y2L_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4Y2R;
   UInt_t          DST_G4HIT_H4Y2R_fUniqueID;
   UInt_t          DST_G4HIT_H4Y2R_fBits;
   Int_t           DST_G4HIT_H4Y2R_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4Y2R_hitmap;
   set<unsigned int> DST_G4HIT_H4Y2R_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4Y2R_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4B;
   UInt_t          DST_G4HIT_H4B_fUniqueID;
   UInt_t          DST_G4HIT_H4B_fBits;
   Int_t           DST_G4HIT_H4B_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4B_hitmap;
   set<unsigned int> DST_G4HIT_H4B_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4B_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_H4T;
   UInt_t          DST_G4HIT_H4T_fUniqueID;
   UInt_t          DST_G4HIT_H4T_fBits;
   Int_t           DST_G4HIT_H4T_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_H4T_hitmap;
   set<unsigned int> DST_G4HIT_H4T_layers;
   map<unsigned int,unsigned int> DST_G4HIT_H4T_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P1Y1;
   UInt_t          DST_G4HIT_P1Y1_fUniqueID;
   UInt_t          DST_G4HIT_P1Y1_fBits;
   Int_t           DST_G4HIT_P1Y1_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P1Y1_hitmap;
   set<unsigned int> DST_G4HIT_P1Y1_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P1Y1_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P1Y2;
   UInt_t          DST_G4HIT_P1Y2_fUniqueID;
   UInt_t          DST_G4HIT_P1Y2_fBits;
   Int_t           DST_G4HIT_P1Y2_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P1Y2_hitmap;
   set<unsigned int> DST_G4HIT_P1Y2_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P1Y2_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P1X1;
   UInt_t          DST_G4HIT_P1X1_fUniqueID;
   UInt_t          DST_G4HIT_P1X1_fBits;
   Int_t           DST_G4HIT_P1X1_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P1X1_hitmap;
   set<unsigned int> DST_G4HIT_P1X1_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P1X1_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P1X2;
   UInt_t          DST_G4HIT_P1X2_fUniqueID;
   UInt_t          DST_G4HIT_P1X2_fBits;
   Int_t           DST_G4HIT_P1X2_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P1X2_hitmap;
   set<unsigned int> DST_G4HIT_P1X2_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P1X2_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P2X1;
   UInt_t          DST_G4HIT_P2X1_fUniqueID;
   UInt_t          DST_G4HIT_P2X1_fBits;
   Int_t           DST_G4HIT_P2X1_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P2X1_hitmap;
   set<unsigned int> DST_G4HIT_P2X1_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P2X1_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P2X2;
   UInt_t          DST_G4HIT_P2X2_fUniqueID;
   UInt_t          DST_G4HIT_P2X2_fBits;
   Int_t           DST_G4HIT_P2X2_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P2X2_hitmap;
   set<unsigned int> DST_G4HIT_P2X2_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P2X2_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P2Y1;
   UInt_t          DST_G4HIT_P2Y1_fUniqueID;
   UInt_t          DST_G4HIT_P2Y1_fBits;
   Int_t           DST_G4HIT_P2Y1_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P2Y1_hitmap;
   set<unsigned int> DST_G4HIT_P2Y1_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P2Y1_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_P2Y2;
   UInt_t          DST_G4HIT_P2Y2_fUniqueID;
   UInt_t          DST_G4HIT_P2Y2_fBits;
   Int_t           DST_G4HIT_P2Y2_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_P2Y2_hitmap;
   set<unsigned int> DST_G4HIT_P2Y2_layers;
   map<unsigned int,unsigned int> DST_G4HIT_P2Y2_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP1TL;
   UInt_t          DST_G4HIT_DP1TL_fUniqueID;
   UInt_t          DST_G4HIT_DP1TL_fBits;
   Int_t           DST_G4HIT_DP1TL_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP1TL_hitmap;
   set<unsigned int> DST_G4HIT_DP1TL_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP1TL_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP1TR;
   UInt_t          DST_G4HIT_DP1TR_fUniqueID;
   UInt_t          DST_G4HIT_DP1TR_fBits;
   Int_t           DST_G4HIT_DP1TR_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP1TR_hitmap;
   set<unsigned int> DST_G4HIT_DP1TR_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP1TR_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP1BL;
   UInt_t          DST_G4HIT_DP1BL_fUniqueID;
   UInt_t          DST_G4HIT_DP1BL_fBits;
   Int_t           DST_G4HIT_DP1BL_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP1BL_hitmap;
   set<unsigned int> DST_G4HIT_DP1BL_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP1BL_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP1BR;
   UInt_t          DST_G4HIT_DP1BR_fUniqueID;
   UInt_t          DST_G4HIT_DP1BR_fBits;
   Int_t           DST_G4HIT_DP1BR_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP1BR_hitmap;
   set<unsigned int> DST_G4HIT_DP1BR_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP1BR_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP2TL;
   UInt_t          DST_G4HIT_DP2TL_fUniqueID;
   UInt_t          DST_G4HIT_DP2TL_fBits;
   Int_t           DST_G4HIT_DP2TL_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP2TL_hitmap;
   set<unsigned int> DST_G4HIT_DP2TL_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP2TL_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP2TR;
   UInt_t          DST_G4HIT_DP2TR_fUniqueID;
   UInt_t          DST_G4HIT_DP2TR_fBits;
   Int_t           DST_G4HIT_DP2TR_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP2TR_hitmap;
   set<unsigned int> DST_G4HIT_DP2TR_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP2TR_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP2BL;
   UInt_t          DST_G4HIT_DP2BL_fUniqueID;
   UInt_t          DST_G4HIT_DP2BL_fBits;
   Int_t           DST_G4HIT_DP2BL_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP2BL_hitmap;
   set<unsigned int> DST_G4HIT_DP2BL_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP2BL_layerMaxID;
 //PHG4HitContainer *DST_G4HIT_DP2BR;
   UInt_t          DST_G4HIT_DP2BR_fUniqueID;
   UInt_t          DST_G4HIT_DP2BR_fBits;
   Int_t           DST_G4HIT_DP2BR_id;
   map<unsigned int,PHG4Hit*> DST_G4HIT_DP2BR_hitmap;
   set<unsigned int> DST_G4HIT_DP2BR_layers;
   map<unsigned int,unsigned int> DST_G4HIT_DP2BR_layerMaxID;
 //PHG4TruthInfoContainer *DST_G4TruthInfo;
   UInt_t          DST_G4TruthInfo_fUniqueID;
   UInt_t          DST_G4TruthInfo_fBits;
   map<int,PHG4Particle*> DST_G4TruthInfo_particlemap;
   map<int,PHG4VtxPoint*> DST_G4TruthInfo_vtxmap;
   map<int,PHG4Shower*> DST_G4TruthInfo_showermap;
   map<int,int>    DST_G4TruthInfo_particle_embed_flags;
   map<int,int>    DST_G4TruthInfo_vertex_embed_flags;
 //SQHitVector_v1  *DST_SQHitVector;
   UInt_t          DST_SQHitVector_fUniqueID;
   UInt_t          DST_SQHitVector_fBits;
   vector<SQHit*>  DST_SQHitVector__vector;
 //SQEvent_v1      *DST_SQEvent;
   UInt_t          DST_SQEvent_fUniqueID;
   UInt_t          DST_SQEvent_fBits;
   Int_t           DST_SQEvent__run_id;
   Int_t           DST_SQEvent__spill_id;
   Int_t           DST_SQEvent__event_id;
   Int_t           DST_SQEvent__coda_event_id;
   UShort_t        DST_SQEvent__trigger;
   Int_t           DST_SQEvent__raw_matrix[5];
   Int_t           DST_SQEvent__after_inh_matrix[5];
   Int_t           DST_SQEvent__data_quality;
   Int_t           DST_SQEvent__vme_time;
   Int_t           DST_SQEvent__qie_presums[4];
   Int_t           DST_SQEvent__qie_trig_cnt;
   Int_t           DST_SQEvent__qie_turn_id;
   Int_t           DST_SQEvent__qie_rf_id;
   Int_t           DST_SQEvent__qie_rf_inte[33];
   Short_t         DST_SQEvent__flag_v1495;
   Short_t         DST_SQEvent__n_board_qie;
   Short_t         DST_SQEvent__n_board_v1495;
   Short_t         DST_SQEvent__n_board_taiwan;
   Short_t         DST_SQEvent__n_board_trig_b;
   Short_t         DST_SQEvent__n_board_trig_c;
 //SQMCEvent_v1    *DST_SQMCEvent;
   UInt_t          DST_SQMCEvent_fUniqueID;
   UInt_t          DST_SQMCEvent_fBits;
   Int_t           DST_SQMCEvent__proc_id;
   Double_t        DST_SQMCEvent__xsec;
   Double_t        DST_SQMCEvent__weight;
   Int_t           DST_SQMCEvent__par_id[4];
   TLorentzVector  DST_SQMCEvent__par_mom[4];
 //SQTrackVector_v1 *DST_SQTruthTrackVector;
   UInt_t          DST_SQTruthTrackVector_fUniqueID;
   UInt_t          DST_SQTruthTrackVector_fBits;
   vector<SQTrack*> DST_SQTruthTrackVector__vector;
 //SQDimuonVector_v1 *DST_SQTruthDimuonVector;
   UInt_t          DST_SQTruthDimuonVector_fUniqueID;
   UInt_t          DST_SQTruthDimuonVector_fBits;
   vector<SQDimuon*> DST_SQTruthDimuonVector__vector;
 //SRecEvent       *DST_SRecEvent;
   UInt_t          DST_SRecEvent_fUniqueID;
   UInt_t          DST_SRecEvent_fBits;
   Short_t         DST_SRecEvent_fRecStatus;
   Int_t           DST_SRecEvent_fRunID;
   Int_t           DST_SRecEvent_fSpillID;
   Int_t           DST_SRecEvent_fEventID;
   Int_t           DST_SRecEvent_fTargetPos;
   Int_t           DST_SRecEvent_fTriggerBits;
   Int_t           DST_SRecEvent_fAllTracks_;
   UInt_t          DST_SRecEvent_fAllTracks_fUniqueID[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   UInt_t          DST_SRecEvent_fAllTracks_fBits[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fChisq[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   vector<int>     DST_SRecEvent_fAllTracks_fHitIndex[kMaxDST_SRecEvent_fAllTracks];
   vector<TMatrixT<double> > DST_SRecEvent_fAllTracks_fState[kMaxDST_SRecEvent_fAllTracks];
   vector<TMatrixT<double> > DST_SRecEvent_fAllTracks_fCovar[kMaxDST_SRecEvent_fAllTracks];
   vector<double>  DST_SRecEvent_fAllTracks_fZ[kMaxDST_SRecEvent_fAllTracks];
   vector<double>  DST_SRecEvent_fAllTracks_fChisqAtNode[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fDumpFacePos[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fDumpPos[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fTargetPos[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fXVertexPos[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fYVertexPos[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fDumpFaceMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fDumpMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fTargetMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fXVertexMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fYVertexMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fVertexMom[kMaxDST_SRecEvent_fAllTracks];
   TVector3        DST_SRecEvent_fAllTracks_fVertexPos[kMaxDST_SRecEvent_fAllTracks];
   Double_t        DST_SRecEvent_fAllTracks_fChisqVertex[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   TMatrixT<double> DST_SRecEvent_fAllTracks_fStateVertex[kMaxDST_SRecEvent_fAllTracks];
   TMatrixT<double> DST_SRecEvent_fAllTracks_fCovarVertex[kMaxDST_SRecEvent_fAllTracks];
   Int_t           DST_SRecEvent_fAllTracks_fKalmanStatus[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Int_t           DST_SRecEvent_fAllTracks_fTriggerID[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Int_t           DST_SRecEvent_fAllTracks_fNPropHitsX[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Int_t           DST_SRecEvent_fAllTracks_fNPropHitsY[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fPropSlopeX[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fPropSlopeY[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fChisqTarget[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fChisqDump[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   Double_t        DST_SRecEvent_fAllTracks_fChisqUpstream[kMaxDST_SRecEvent_fAllTracks];   //[DST.SRecEvent.fAllTracks_]
   vector<TVector3> DST_SRecEvent_fAllTracks_fGFDetPlaneVec[3][kMaxDST_SRecEvent_fAllTracks];
   vector<TVectorT<double> > DST_SRecEvent_fAllTracks_fGFAuxInfo[kMaxDST_SRecEvent_fAllTracks];
   vector<TVectorT<double> > DST_SRecEvent_fAllTracks_fGFStateVec[kMaxDST_SRecEvent_fAllTracks];
   vector<TMatrixTSym<double> > DST_SRecEvent_fAllTracks_fGFCov[kMaxDST_SRecEvent_fAllTracks];
   Int_t           DST_SRecEvent_fDimuons_;
   UInt_t          DST_SRecEvent_fDimuons_fUniqueID[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   UInt_t          DST_SRecEvent_fDimuons_fBits[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Int_t           DST_SRecEvent_fDimuons_trackID_pos[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Int_t           DST_SRecEvent_fDimuons_trackID_neg[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   TLorentzVector  DST_SRecEvent_fDimuons_p_pos[kMaxDST_SRecEvent_fDimuons];
   TLorentzVector  DST_SRecEvent_fDimuons_p_neg[kMaxDST_SRecEvent_fDimuons];
   TLorentzVector  DST_SRecEvent_fDimuons_p_pos_single[kMaxDST_SRecEvent_fDimuons];
   TLorentzVector  DST_SRecEvent_fDimuons_p_neg_single[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_vtx[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_vtx_pos[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_vtx_neg[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_proj_target_pos[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_proj_dump_pos[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_proj_target_neg[kMaxDST_SRecEvent_fDimuons];
   TVector3        DST_SRecEvent_fDimuons_proj_dump_neg[kMaxDST_SRecEvent_fDimuons];
   Double_t        DST_SRecEvent_fDimuons_mass[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_pT[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_xF[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_x1[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_x2[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_costh[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_phi[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_mass_single[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_single[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_kf[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_vx[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_target[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_dump[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   Double_t        DST_SRecEvent_fDimuons_chisq_upstream[kMaxDST_SRecEvent_fDimuons];   //[DST.SRecEvent.fDimuons_]
   map<int,int>    DST_SRecEvent_fLocalID;
   Int_t           DST_SRecEvent_fSource1;
   Int_t           DST_SRecEvent_fSource2;

   // List of branches
   TBranch        *b_DST_PHHepMCGenEventMap_fUniqueID;   //!
   TBranch        *b_DST_PHHepMCGenEventMap_fBits;   //!
   TBranch        *b_DST_PHHepMCGenEventMap__map;   //!
   TBranch        *b_DST_PHG4INEVENT_fUniqueID;   //!
   TBranch        *b_DST_PHG4INEVENT_fBits;   //!
   TBranch        *b_DST_PHG4INEVENT_vtxlist;   //!
   TBranch        *b_DST_PHG4INEVENT_particlelist;   //!
   TBranch        *b_DST_PHG4INEVENT_embedded_particlelist;   //!
   TBranch        *b_DST_G4HIT_D0U_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0U_fBits;   //!
   TBranch        *b_DST_G4HIT_D0U_id;   //!
   TBranch        *b_DST_G4HIT_D0U_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0U_layers;   //!
   TBranch        *b_DST_G4HIT_D0U_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D0Up_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0Up_fBits;   //!
   TBranch        *b_DST_G4HIT_D0Up_id;   //!
   TBranch        *b_DST_G4HIT_D0Up_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0Up_layers;   //!
   TBranch        *b_DST_G4HIT_D0Up_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D0X_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0X_fBits;   //!
   TBranch        *b_DST_G4HIT_D0X_id;   //!
   TBranch        *b_DST_G4HIT_D0X_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0X_layers;   //!
   TBranch        *b_DST_G4HIT_D0X_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D0Xp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0Xp_fBits;   //!
   TBranch        *b_DST_G4HIT_D0Xp_id;   //!
   TBranch        *b_DST_G4HIT_D0Xp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0Xp_layers;   //!
   TBranch        *b_DST_G4HIT_D0Xp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D0V_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0V_fBits;   //!
   TBranch        *b_DST_G4HIT_D0V_id;   //!
   TBranch        *b_DST_G4HIT_D0V_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0V_layers;   //!
   TBranch        *b_DST_G4HIT_D0V_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D0Vp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D0Vp_fBits;   //!
   TBranch        *b_DST_G4HIT_D0Vp_id;   //!
   TBranch        *b_DST_G4HIT_D0Vp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D0Vp_layers;   //!
   TBranch        *b_DST_G4HIT_D0Vp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2V_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2V_fBits;   //!
   TBranch        *b_DST_G4HIT_D2V_id;   //!
   TBranch        *b_DST_G4HIT_D2V_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2V_layers;   //!
   TBranch        *b_DST_G4HIT_D2V_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2Vp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2Vp_fBits;   //!
   TBranch        *b_DST_G4HIT_D2Vp_id;   //!
   TBranch        *b_DST_G4HIT_D2Vp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2Vp_layers;   //!
   TBranch        *b_DST_G4HIT_D2Vp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2Xp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2Xp_fBits;   //!
   TBranch        *b_DST_G4HIT_D2Xp_id;   //!
   TBranch        *b_DST_G4HIT_D2Xp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2Xp_layers;   //!
   TBranch        *b_DST_G4HIT_D2Xp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2X_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2X_fBits;   //!
   TBranch        *b_DST_G4HIT_D2X_id;   //!
   TBranch        *b_DST_G4HIT_D2X_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2X_layers;   //!
   TBranch        *b_DST_G4HIT_D2X_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2U_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2U_fBits;   //!
   TBranch        *b_DST_G4HIT_D2U_id;   //!
   TBranch        *b_DST_G4HIT_D2U_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2U_layers;   //!
   TBranch        *b_DST_G4HIT_D2U_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D2Up_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D2Up_fBits;   //!
   TBranch        *b_DST_G4HIT_D2Up_id;   //!
   TBranch        *b_DST_G4HIT_D2Up_hitmap;   //!
   TBranch        *b_DST_G4HIT_D2Up_layers;   //!
   TBranch        *b_DST_G4HIT_D2Up_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pVp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pVp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pVp_id;   //!
   TBranch        *b_DST_G4HIT_D3pVp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pVp_layers;   //!
   TBranch        *b_DST_G4HIT_D3pVp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pV_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pV_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pV_id;   //!
   TBranch        *b_DST_G4HIT_D3pV_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pV_layers;   //!
   TBranch        *b_DST_G4HIT_D3pV_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pXp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pXp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pXp_id;   //!
   TBranch        *b_DST_G4HIT_D3pXp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pXp_layers;   //!
   TBranch        *b_DST_G4HIT_D3pXp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pX_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pX_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pX_id;   //!
   TBranch        *b_DST_G4HIT_D3pX_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pX_layers;   //!
   TBranch        *b_DST_G4HIT_D3pX_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pUp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pUp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pUp_id;   //!
   TBranch        *b_DST_G4HIT_D3pUp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pUp_layers;   //!
   TBranch        *b_DST_G4HIT_D3pUp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3pU_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3pU_fBits;   //!
   TBranch        *b_DST_G4HIT_D3pU_id;   //!
   TBranch        *b_DST_G4HIT_D3pU_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3pU_layers;   //!
   TBranch        *b_DST_G4HIT_D3pU_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mVp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mVp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mVp_id;   //!
   TBranch        *b_DST_G4HIT_D3mVp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mVp_layers;   //!
   TBranch        *b_DST_G4HIT_D3mVp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mV_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mV_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mV_id;   //!
   TBranch        *b_DST_G4HIT_D3mV_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mV_layers;   //!
   TBranch        *b_DST_G4HIT_D3mV_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mXp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mXp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mXp_id;   //!
   TBranch        *b_DST_G4HIT_D3mXp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mXp_layers;   //!
   TBranch        *b_DST_G4HIT_D3mXp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mX_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mX_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mX_id;   //!
   TBranch        *b_DST_G4HIT_D3mX_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mX_layers;   //!
   TBranch        *b_DST_G4HIT_D3mX_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mUp_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mUp_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mUp_id;   //!
   TBranch        *b_DST_G4HIT_D3mUp_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mUp_layers;   //!
   TBranch        *b_DST_G4HIT_D3mUp_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_D3mU_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_D3mU_fBits;   //!
   TBranch        *b_DST_G4HIT_D3mU_id;   //!
   TBranch        *b_DST_G4HIT_D3mU_hitmap;   //!
   TBranch        *b_DST_G4HIT_D3mU_layers;   //!
   TBranch        *b_DST_G4HIT_D3mU_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H1B_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H1B_fBits;   //!
   TBranch        *b_DST_G4HIT_H1B_id;   //!
   TBranch        *b_DST_G4HIT_H1B_hitmap;   //!
   TBranch        *b_DST_G4HIT_H1B_layers;   //!
   TBranch        *b_DST_G4HIT_H1B_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H1T_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H1T_fBits;   //!
   TBranch        *b_DST_G4HIT_H1T_id;   //!
   TBranch        *b_DST_G4HIT_H1T_hitmap;   //!
   TBranch        *b_DST_G4HIT_H1T_layers;   //!
   TBranch        *b_DST_G4HIT_H1T_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H1L_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H1L_fBits;   //!
   TBranch        *b_DST_G4HIT_H1L_id;   //!
   TBranch        *b_DST_G4HIT_H1L_hitmap;   //!
   TBranch        *b_DST_G4HIT_H1L_layers;   //!
   TBranch        *b_DST_G4HIT_H1L_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H1R_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H1R_fBits;   //!
   TBranch        *b_DST_G4HIT_H1R_id;   //!
   TBranch        *b_DST_G4HIT_H1R_hitmap;   //!
   TBranch        *b_DST_G4HIT_H1R_layers;   //!
   TBranch        *b_DST_G4HIT_H1R_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H2L_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H2L_fBits;   //!
   TBranch        *b_DST_G4HIT_H2L_id;   //!
   TBranch        *b_DST_G4HIT_H2L_hitmap;   //!
   TBranch        *b_DST_G4HIT_H2L_layers;   //!
   TBranch        *b_DST_G4HIT_H2L_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H2R_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H2R_fBits;   //!
   TBranch        *b_DST_G4HIT_H2R_id;   //!
   TBranch        *b_DST_G4HIT_H2R_hitmap;   //!
   TBranch        *b_DST_G4HIT_H2R_layers;   //!
   TBranch        *b_DST_G4HIT_H2R_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H2B_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H2B_fBits;   //!
   TBranch        *b_DST_G4HIT_H2B_id;   //!
   TBranch        *b_DST_G4HIT_H2B_hitmap;   //!
   TBranch        *b_DST_G4HIT_H2B_layers;   //!
   TBranch        *b_DST_G4HIT_H2B_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H2T_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H2T_fBits;   //!
   TBranch        *b_DST_G4HIT_H2T_id;   //!
   TBranch        *b_DST_G4HIT_H2T_hitmap;   //!
   TBranch        *b_DST_G4HIT_H2T_layers;   //!
   TBranch        *b_DST_G4HIT_H2T_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H3B_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H3B_fBits;   //!
   TBranch        *b_DST_G4HIT_H3B_id;   //!
   TBranch        *b_DST_G4HIT_H3B_hitmap;   //!
   TBranch        *b_DST_G4HIT_H3B_layers;   //!
   TBranch        *b_DST_G4HIT_H3B_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H3T_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H3T_fBits;   //!
   TBranch        *b_DST_G4HIT_H3T_id;   //!
   TBranch        *b_DST_G4HIT_H3T_hitmap;   //!
   TBranch        *b_DST_G4HIT_H3T_layers;   //!
   TBranch        *b_DST_G4HIT_H3T_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_fBits;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_id;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_layers;   //!
   TBranch        *b_DST_G4HIT_H4Y1L_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_fBits;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_id;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_layers;   //!
   TBranch        *b_DST_G4HIT_H4Y1R_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_fBits;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_id;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_layers;   //!
   TBranch        *b_DST_G4HIT_H4Y2L_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_fBits;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_id;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_layers;   //!
   TBranch        *b_DST_G4HIT_H4Y2R_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4B_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4B_fBits;   //!
   TBranch        *b_DST_G4HIT_H4B_id;   //!
   TBranch        *b_DST_G4HIT_H4B_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4B_layers;   //!
   TBranch        *b_DST_G4HIT_H4B_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_H4T_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_H4T_fBits;   //!
   TBranch        *b_DST_G4HIT_H4T_id;   //!
   TBranch        *b_DST_G4HIT_H4T_hitmap;   //!
   TBranch        *b_DST_G4HIT_H4T_layers;   //!
   TBranch        *b_DST_G4HIT_H4T_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P1Y1_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P1Y1_fBits;   //!
   TBranch        *b_DST_G4HIT_P1Y1_id;   //!
   TBranch        *b_DST_G4HIT_P1Y1_hitmap;   //!
   TBranch        *b_DST_G4HIT_P1Y1_layers;   //!
   TBranch        *b_DST_G4HIT_P1Y1_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P1Y2_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P1Y2_fBits;   //!
   TBranch        *b_DST_G4HIT_P1Y2_id;   //!
   TBranch        *b_DST_G4HIT_P1Y2_hitmap;   //!
   TBranch        *b_DST_G4HIT_P1Y2_layers;   //!
   TBranch        *b_DST_G4HIT_P1Y2_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P1X1_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P1X1_fBits;   //!
   TBranch        *b_DST_G4HIT_P1X1_id;   //!
   TBranch        *b_DST_G4HIT_P1X1_hitmap;   //!
   TBranch        *b_DST_G4HIT_P1X1_layers;   //!
   TBranch        *b_DST_G4HIT_P1X1_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P1X2_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P1X2_fBits;   //!
   TBranch        *b_DST_G4HIT_P1X2_id;   //!
   TBranch        *b_DST_G4HIT_P1X2_hitmap;   //!
   TBranch        *b_DST_G4HIT_P1X2_layers;   //!
   TBranch        *b_DST_G4HIT_P1X2_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P2X1_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P2X1_fBits;   //!
   TBranch        *b_DST_G4HIT_P2X1_id;   //!
   TBranch        *b_DST_G4HIT_P2X1_hitmap;   //!
   TBranch        *b_DST_G4HIT_P2X1_layers;   //!
   TBranch        *b_DST_G4HIT_P2X1_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P2X2_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P2X2_fBits;   //!
   TBranch        *b_DST_G4HIT_P2X2_id;   //!
   TBranch        *b_DST_G4HIT_P2X2_hitmap;   //!
   TBranch        *b_DST_G4HIT_P2X2_layers;   //!
   TBranch        *b_DST_G4HIT_P2X2_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P2Y1_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P2Y1_fBits;   //!
   TBranch        *b_DST_G4HIT_P2Y1_id;   //!
   TBranch        *b_DST_G4HIT_P2Y1_hitmap;   //!
   TBranch        *b_DST_G4HIT_P2Y1_layers;   //!
   TBranch        *b_DST_G4HIT_P2Y1_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_P2Y2_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_P2Y2_fBits;   //!
   TBranch        *b_DST_G4HIT_P2Y2_id;   //!
   TBranch        *b_DST_G4HIT_P2Y2_hitmap;   //!
   TBranch        *b_DST_G4HIT_P2Y2_layers;   //!
   TBranch        *b_DST_G4HIT_P2Y2_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP1TL_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP1TL_fBits;   //!
   TBranch        *b_DST_G4HIT_DP1TL_id;   //!
   TBranch        *b_DST_G4HIT_DP1TL_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP1TL_layers;   //!
   TBranch        *b_DST_G4HIT_DP1TL_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP1TR_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP1TR_fBits;   //!
   TBranch        *b_DST_G4HIT_DP1TR_id;   //!
   TBranch        *b_DST_G4HIT_DP1TR_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP1TR_layers;   //!
   TBranch        *b_DST_G4HIT_DP1TR_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP1BL_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP1BL_fBits;   //!
   TBranch        *b_DST_G4HIT_DP1BL_id;   //!
   TBranch        *b_DST_G4HIT_DP1BL_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP1BL_layers;   //!
   TBranch        *b_DST_G4HIT_DP1BL_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP1BR_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP1BR_fBits;   //!
   TBranch        *b_DST_G4HIT_DP1BR_id;   //!
   TBranch        *b_DST_G4HIT_DP1BR_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP1BR_layers;   //!
   TBranch        *b_DST_G4HIT_DP1BR_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP2TL_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP2TL_fBits;   //!
   TBranch        *b_DST_G4HIT_DP2TL_id;   //!
   TBranch        *b_DST_G4HIT_DP2TL_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP2TL_layers;   //!
   TBranch        *b_DST_G4HIT_DP2TL_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP2TR_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP2TR_fBits;   //!
   TBranch        *b_DST_G4HIT_DP2TR_id;   //!
   TBranch        *b_DST_G4HIT_DP2TR_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP2TR_layers;   //!
   TBranch        *b_DST_G4HIT_DP2TR_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP2BL_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP2BL_fBits;   //!
   TBranch        *b_DST_G4HIT_DP2BL_id;   //!
   TBranch        *b_DST_G4HIT_DP2BL_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP2BL_layers;   //!
   TBranch        *b_DST_G4HIT_DP2BL_layerMaxID;   //!
   TBranch        *b_DST_G4HIT_DP2BR_fUniqueID;   //!
   TBranch        *b_DST_G4HIT_DP2BR_fBits;   //!
   TBranch        *b_DST_G4HIT_DP2BR_id;   //!
   TBranch        *b_DST_G4HIT_DP2BR_hitmap;   //!
   TBranch        *b_DST_G4HIT_DP2BR_layers;   //!
   TBranch        *b_DST_G4HIT_DP2BR_layerMaxID;   //!
   TBranch        *b_DST_G4TruthInfo_fUniqueID;   //!
   TBranch        *b_DST_G4TruthInfo_fBits;   //!
   TBranch        *b_DST_G4TruthInfo_particlemap;   //!
   TBranch        *b_DST_G4TruthInfo_vtxmap;   //!
   TBranch        *b_DST_G4TruthInfo_showermap;   //!
   TBranch        *b_DST_G4TruthInfo_particle_embed_flags;   //!
   TBranch        *b_DST_G4TruthInfo_vertex_embed_flags;   //!
   TBranch        *b_DST_SQHitVector_fUniqueID;   //!
   TBranch        *b_DST_SQHitVector_fBits;   //!
   TBranch        *b_DST_SQHitVector__vector;   //!
   TBranch        *b_DST_SQEvent_fUniqueID;   //!
   TBranch        *b_DST_SQEvent_fBits;   //!
   TBranch        *b_DST_SQEvent__run_id;   //!
   TBranch        *b_DST_SQEvent__spill_id;   //!
   TBranch        *b_DST_SQEvent__event_id;   //!
   TBranch        *b_DST_SQEvent__coda_event_id;   //!
   TBranch        *b_DST_SQEvent__trigger;   //!
   TBranch        *b_DST_SQEvent__raw_matrix;   //!
   TBranch        *b_DST_SQEvent__after_inh_matrix;   //!
   TBranch        *b_DST_SQEvent__data_quality;   //!
   TBranch        *b_DST_SQEvent__vme_time;   //!
   TBranch        *b_DST_SQEvent__qie_presums;   //!
   TBranch        *b_DST_SQEvent__qie_trig_cnt;   //!
   TBranch        *b_DST_SQEvent__qie_turn_id;   //!
   TBranch        *b_DST_SQEvent__qie_rf_id;   //!
   TBranch        *b_DST_SQEvent__qie_rf_inte;   //!
   TBranch        *b_DST_SQEvent__flag_v1495;   //!
   TBranch        *b_DST_SQEvent__n_board_qie;   //!
   TBranch        *b_DST_SQEvent__n_board_v1495;   //!
   TBranch        *b_DST_SQEvent__n_board_taiwan;   //!
   TBranch        *b_DST_SQEvent__n_board_trig_b;   //!
   TBranch        *b_DST_SQEvent__n_board_trig_c;   //!
   TBranch        *b_DST_SQMCEvent_fUniqueID;   //!
   TBranch        *b_DST_SQMCEvent_fBits;   //!
   TBranch        *b_DST_SQMCEvent__proc_id;   //!
   TBranch        *b_DST_SQMCEvent__xsec;   //!
   TBranch        *b_DST_SQMCEvent__weight;   //!
   TBranch        *b_DST_SQMCEvent__par_id;   //!
   TBranch        *b_DST_SQMCEvent__par_mom;   //!
   TBranch        *b_DST_SQTruthTrackVector_fUniqueID;   //!
   TBranch        *b_DST_SQTruthTrackVector_fBits;   //!
   TBranch        *b_DST_SQTruthTrackVector__vector;   //!
   TBranch        *b_DST_SQTruthDimuonVector_fUniqueID;   //!
   TBranch        *b_DST_SQTruthDimuonVector_fBits;   //!
   TBranch        *b_DST_SQTruthDimuonVector__vector;   //!
   TBranch        *b_DST_SRecEvent_fUniqueID;   //!
   TBranch        *b_DST_SRecEvent_fBits;   //!
   TBranch        *b_DST_SRecEvent_fRecStatus;   //!
   TBranch        *b_DST_SRecEvent_fRunID;   //!
   TBranch        *b_DST_SRecEvent_fSpillID;   //!
   TBranch        *b_DST_SRecEvent_fEventID;   //!
   TBranch        *b_DST_SRecEvent_fTargetPos;   //!
   TBranch        *b_DST_SRecEvent_fTriggerBits;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fUniqueID;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fBits;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisq;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fHitIndex;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fState;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fCovar;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fZ;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisqAtNode;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fDumpFacePos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fDumpPos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fTargetPos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fXVertexPos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fYVertexPos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fDumpFaceMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fDumpMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fTargetMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fXVertexMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fYVertexMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fVertexMom;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fVertexPos;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisqVertex;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fStateVertex;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fCovarVertex;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fKalmanStatus;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fTriggerID;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fNPropHitsX;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fNPropHitsY;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fPropSlopeX;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fPropSlopeY;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisqTarget;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisqDump;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fChisqUpstream;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fGFDetPlaneVec;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fGFAuxInfo;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fGFStateVec;   //!
   TBranch        *b_DST_SRecEvent_fAllTracks_fGFCov;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_fUniqueID;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_fBits;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_trackID_pos;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_trackID_neg;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_p_pos;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_p_neg;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_p_pos_single;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_p_neg_single;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_vtx;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_vtx_pos;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_vtx_neg;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_proj_target_pos;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_proj_dump_pos;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_proj_target_neg;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_proj_dump_neg;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_mass;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_pT;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_xF;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_x1;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_x2;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_costh;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_phi;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_mass_single;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_single;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_kf;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_vx;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_target;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_dump;   //!
   TBranch        *b_DST_SRecEvent_fDimuons_chisq_upstream;   //!
   TBranch        *b_DST_SRecEvent_fLocalID;   //!
   TBranch        *b_DST_SRecEvent_fSource1;   //!
   TBranch        *b_DST_SRecEvent_fSource2;   //!

   e1039data_tree(TTree *tree=0);
   virtual ~e1039data_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef e1039data_tree_cxx
e1039data_tree::e1039data_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DST.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DST.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

e1039data_tree::~e1039data_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t e1039data_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t e1039data_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void e1039data_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("DST.PHHepMCGenEventMap.fUniqueID", &DST_PHHepMCGenEventMap_fUniqueID, &b_DST_PHHepMCGenEventMap_fUniqueID);
   fChain->SetBranchAddress("DST.PHHepMCGenEventMap.fBits", &DST_PHHepMCGenEventMap_fBits, &b_DST_PHHepMCGenEventMap_fBits);
   fChain->SetBranchAddress("DST.PHHepMCGenEventMap._map", &DST_PHHepMCGenEventMap__map, &b_DST_PHHepMCGenEventMap__map);
   fChain->SetBranchAddress("DST.PHG4INEVENT.fUniqueID", &DST_PHG4INEVENT_fUniqueID, &b_DST_PHG4INEVENT_fUniqueID);
   fChain->SetBranchAddress("DST.PHG4INEVENT.fBits", &DST_PHG4INEVENT_fBits, &b_DST_PHG4INEVENT_fBits);
   fChain->SetBranchAddress("DST.PHG4INEVENT.vtxlist", &DST_PHG4INEVENT_vtxlist, &b_DST_PHG4INEVENT_vtxlist);
   fChain->SetBranchAddress("DST.PHG4INEVENT.particlelist", &DST_PHG4INEVENT_particlelist, &b_DST_PHG4INEVENT_particlelist);
   fChain->SetBranchAddress("DST.PHG4INEVENT.embedded_particlelist", &DST_PHG4INEVENT_embedded_particlelist, &b_DST_PHG4INEVENT_embedded_particlelist);
   fChain->SetBranchAddress("DST.G4HIT_D0U.fUniqueID", &DST_G4HIT_D0U_fUniqueID, &b_DST_G4HIT_D0U_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0U.fBits", &DST_G4HIT_D0U_fBits, &b_DST_G4HIT_D0U_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0U.id", &DST_G4HIT_D0U_id, &b_DST_G4HIT_D0U_id);
   fChain->SetBranchAddress("DST.G4HIT_D0U.hitmap", &DST_G4HIT_D0U_hitmap, &b_DST_G4HIT_D0U_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0U.layers", &DST_G4HIT_D0U_layers, &b_DST_G4HIT_D0U_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0U.layerMaxID", &DST_G4HIT_D0U_layerMaxID, &b_DST_G4HIT_D0U_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.fUniqueID", &DST_G4HIT_D0Up_fUniqueID, &b_DST_G4HIT_D0Up_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.fBits", &DST_G4HIT_D0Up_fBits, &b_DST_G4HIT_D0Up_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.id", &DST_G4HIT_D0Up_id, &b_DST_G4HIT_D0Up_id);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.hitmap", &DST_G4HIT_D0Up_hitmap, &b_DST_G4HIT_D0Up_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.layers", &DST_G4HIT_D0Up_layers, &b_DST_G4HIT_D0Up_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0Up.layerMaxID", &DST_G4HIT_D0Up_layerMaxID, &b_DST_G4HIT_D0Up_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D0X.fUniqueID", &DST_G4HIT_D0X_fUniqueID, &b_DST_G4HIT_D0X_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0X.fBits", &DST_G4HIT_D0X_fBits, &b_DST_G4HIT_D0X_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0X.id", &DST_G4HIT_D0X_id, &b_DST_G4HIT_D0X_id);
   fChain->SetBranchAddress("DST.G4HIT_D0X.hitmap", &DST_G4HIT_D0X_hitmap, &b_DST_G4HIT_D0X_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0X.layers", &DST_G4HIT_D0X_layers, &b_DST_G4HIT_D0X_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0X.layerMaxID", &DST_G4HIT_D0X_layerMaxID, &b_DST_G4HIT_D0X_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.fUniqueID", &DST_G4HIT_D0Xp_fUniqueID, &b_DST_G4HIT_D0Xp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.fBits", &DST_G4HIT_D0Xp_fBits, &b_DST_G4HIT_D0Xp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.id", &DST_G4HIT_D0Xp_id, &b_DST_G4HIT_D0Xp_id);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.hitmap", &DST_G4HIT_D0Xp_hitmap, &b_DST_G4HIT_D0Xp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.layers", &DST_G4HIT_D0Xp_layers, &b_DST_G4HIT_D0Xp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0Xp.layerMaxID", &DST_G4HIT_D0Xp_layerMaxID, &b_DST_G4HIT_D0Xp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D0V.fUniqueID", &DST_G4HIT_D0V_fUniqueID, &b_DST_G4HIT_D0V_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0V.fBits", &DST_G4HIT_D0V_fBits, &b_DST_G4HIT_D0V_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0V.id", &DST_G4HIT_D0V_id, &b_DST_G4HIT_D0V_id);
   fChain->SetBranchAddress("DST.G4HIT_D0V.hitmap", &DST_G4HIT_D0V_hitmap, &b_DST_G4HIT_D0V_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0V.layers", &DST_G4HIT_D0V_layers, &b_DST_G4HIT_D0V_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0V.layerMaxID", &DST_G4HIT_D0V_layerMaxID, &b_DST_G4HIT_D0V_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.fUniqueID", &DST_G4HIT_D0Vp_fUniqueID, &b_DST_G4HIT_D0Vp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.fBits", &DST_G4HIT_D0Vp_fBits, &b_DST_G4HIT_D0Vp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.id", &DST_G4HIT_D0Vp_id, &b_DST_G4HIT_D0Vp_id);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.hitmap", &DST_G4HIT_D0Vp_hitmap, &b_DST_G4HIT_D0Vp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.layers", &DST_G4HIT_D0Vp_layers, &b_DST_G4HIT_D0Vp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D0Vp.layerMaxID", &DST_G4HIT_D0Vp_layerMaxID, &b_DST_G4HIT_D0Vp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2V.fUniqueID", &DST_G4HIT_D2V_fUniqueID, &b_DST_G4HIT_D2V_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2V.fBits", &DST_G4HIT_D2V_fBits, &b_DST_G4HIT_D2V_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2V.id", &DST_G4HIT_D2V_id, &b_DST_G4HIT_D2V_id);
   fChain->SetBranchAddress("DST.G4HIT_D2V.hitmap", &DST_G4HIT_D2V_hitmap, &b_DST_G4HIT_D2V_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2V.layers", &DST_G4HIT_D2V_layers, &b_DST_G4HIT_D2V_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2V.layerMaxID", &DST_G4HIT_D2V_layerMaxID, &b_DST_G4HIT_D2V_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.fUniqueID", &DST_G4HIT_D2Vp_fUniqueID, &b_DST_G4HIT_D2Vp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.fBits", &DST_G4HIT_D2Vp_fBits, &b_DST_G4HIT_D2Vp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.id", &DST_G4HIT_D2Vp_id, &b_DST_G4HIT_D2Vp_id);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.hitmap", &DST_G4HIT_D2Vp_hitmap, &b_DST_G4HIT_D2Vp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.layers", &DST_G4HIT_D2Vp_layers, &b_DST_G4HIT_D2Vp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2Vp.layerMaxID", &DST_G4HIT_D2Vp_layerMaxID, &b_DST_G4HIT_D2Vp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.fUniqueID", &DST_G4HIT_D2Xp_fUniqueID, &b_DST_G4HIT_D2Xp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.fBits", &DST_G4HIT_D2Xp_fBits, &b_DST_G4HIT_D2Xp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.id", &DST_G4HIT_D2Xp_id, &b_DST_G4HIT_D2Xp_id);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.hitmap", &DST_G4HIT_D2Xp_hitmap, &b_DST_G4HIT_D2Xp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.layers", &DST_G4HIT_D2Xp_layers, &b_DST_G4HIT_D2Xp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2Xp.layerMaxID", &DST_G4HIT_D2Xp_layerMaxID, &b_DST_G4HIT_D2Xp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2X.fUniqueID", &DST_G4HIT_D2X_fUniqueID, &b_DST_G4HIT_D2X_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2X.fBits", &DST_G4HIT_D2X_fBits, &b_DST_G4HIT_D2X_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2X.id", &DST_G4HIT_D2X_id, &b_DST_G4HIT_D2X_id);
   fChain->SetBranchAddress("DST.G4HIT_D2X.hitmap", &DST_G4HIT_D2X_hitmap, &b_DST_G4HIT_D2X_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2X.layers", &DST_G4HIT_D2X_layers, &b_DST_G4HIT_D2X_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2X.layerMaxID", &DST_G4HIT_D2X_layerMaxID, &b_DST_G4HIT_D2X_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2U.fUniqueID", &DST_G4HIT_D2U_fUniqueID, &b_DST_G4HIT_D2U_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2U.fBits", &DST_G4HIT_D2U_fBits, &b_DST_G4HIT_D2U_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2U.id", &DST_G4HIT_D2U_id, &b_DST_G4HIT_D2U_id);
   fChain->SetBranchAddress("DST.G4HIT_D2U.hitmap", &DST_G4HIT_D2U_hitmap, &b_DST_G4HIT_D2U_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2U.layers", &DST_G4HIT_D2U_layers, &b_DST_G4HIT_D2U_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2U.layerMaxID", &DST_G4HIT_D2U_layerMaxID, &b_DST_G4HIT_D2U_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.fUniqueID", &DST_G4HIT_D2Up_fUniqueID, &b_DST_G4HIT_D2Up_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.fBits", &DST_G4HIT_D2Up_fBits, &b_DST_G4HIT_D2Up_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.id", &DST_G4HIT_D2Up_id, &b_DST_G4HIT_D2Up_id);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.hitmap", &DST_G4HIT_D2Up_hitmap, &b_DST_G4HIT_D2Up_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.layers", &DST_G4HIT_D2Up_layers, &b_DST_G4HIT_D2Up_layers);
   fChain->SetBranchAddress("DST.G4HIT_D2Up.layerMaxID", &DST_G4HIT_D2Up_layerMaxID, &b_DST_G4HIT_D2Up_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.fUniqueID", &DST_G4HIT_D3pVp_fUniqueID, &b_DST_G4HIT_D3pVp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.fBits", &DST_G4HIT_D3pVp_fBits, &b_DST_G4HIT_D3pVp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.id", &DST_G4HIT_D3pVp_id, &b_DST_G4HIT_D3pVp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.hitmap", &DST_G4HIT_D3pVp_hitmap, &b_DST_G4HIT_D3pVp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.layers", &DST_G4HIT_D3pVp_layers, &b_DST_G4HIT_D3pVp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pVp.layerMaxID", &DST_G4HIT_D3pVp_layerMaxID, &b_DST_G4HIT_D3pVp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.fUniqueID", &DST_G4HIT_D3pV_fUniqueID, &b_DST_G4HIT_D3pV_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.fBits", &DST_G4HIT_D3pV_fBits, &b_DST_G4HIT_D3pV_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.id", &DST_G4HIT_D3pV_id, &b_DST_G4HIT_D3pV_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.hitmap", &DST_G4HIT_D3pV_hitmap, &b_DST_G4HIT_D3pV_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.layers", &DST_G4HIT_D3pV_layers, &b_DST_G4HIT_D3pV_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pV.layerMaxID", &DST_G4HIT_D3pV_layerMaxID, &b_DST_G4HIT_D3pV_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.fUniqueID", &DST_G4HIT_D3pXp_fUniqueID, &b_DST_G4HIT_D3pXp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.fBits", &DST_G4HIT_D3pXp_fBits, &b_DST_G4HIT_D3pXp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.id", &DST_G4HIT_D3pXp_id, &b_DST_G4HIT_D3pXp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.hitmap", &DST_G4HIT_D3pXp_hitmap, &b_DST_G4HIT_D3pXp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.layers", &DST_G4HIT_D3pXp_layers, &b_DST_G4HIT_D3pXp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pXp.layerMaxID", &DST_G4HIT_D3pXp_layerMaxID, &b_DST_G4HIT_D3pXp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.fUniqueID", &DST_G4HIT_D3pX_fUniqueID, &b_DST_G4HIT_D3pX_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.fBits", &DST_G4HIT_D3pX_fBits, &b_DST_G4HIT_D3pX_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.id", &DST_G4HIT_D3pX_id, &b_DST_G4HIT_D3pX_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.hitmap", &DST_G4HIT_D3pX_hitmap, &b_DST_G4HIT_D3pX_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.layers", &DST_G4HIT_D3pX_layers, &b_DST_G4HIT_D3pX_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pX.layerMaxID", &DST_G4HIT_D3pX_layerMaxID, &b_DST_G4HIT_D3pX_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.fUniqueID", &DST_G4HIT_D3pUp_fUniqueID, &b_DST_G4HIT_D3pUp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.fBits", &DST_G4HIT_D3pUp_fBits, &b_DST_G4HIT_D3pUp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.id", &DST_G4HIT_D3pUp_id, &b_DST_G4HIT_D3pUp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.hitmap", &DST_G4HIT_D3pUp_hitmap, &b_DST_G4HIT_D3pUp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.layers", &DST_G4HIT_D3pUp_layers, &b_DST_G4HIT_D3pUp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pUp.layerMaxID", &DST_G4HIT_D3pUp_layerMaxID, &b_DST_G4HIT_D3pUp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.fUniqueID", &DST_G4HIT_D3pU_fUniqueID, &b_DST_G4HIT_D3pU_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.fBits", &DST_G4HIT_D3pU_fBits, &b_DST_G4HIT_D3pU_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.id", &DST_G4HIT_D3pU_id, &b_DST_G4HIT_D3pU_id);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.hitmap", &DST_G4HIT_D3pU_hitmap, &b_DST_G4HIT_D3pU_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.layers", &DST_G4HIT_D3pU_layers, &b_DST_G4HIT_D3pU_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3pU.layerMaxID", &DST_G4HIT_D3pU_layerMaxID, &b_DST_G4HIT_D3pU_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.fUniqueID", &DST_G4HIT_D3mVp_fUniqueID, &b_DST_G4HIT_D3mVp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.fBits", &DST_G4HIT_D3mVp_fBits, &b_DST_G4HIT_D3mVp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.id", &DST_G4HIT_D3mVp_id, &b_DST_G4HIT_D3mVp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.hitmap", &DST_G4HIT_D3mVp_hitmap, &b_DST_G4HIT_D3mVp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.layers", &DST_G4HIT_D3mVp_layers, &b_DST_G4HIT_D3mVp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mVp.layerMaxID", &DST_G4HIT_D3mVp_layerMaxID, &b_DST_G4HIT_D3mVp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.fUniqueID", &DST_G4HIT_D3mV_fUniqueID, &b_DST_G4HIT_D3mV_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.fBits", &DST_G4HIT_D3mV_fBits, &b_DST_G4HIT_D3mV_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.id", &DST_G4HIT_D3mV_id, &b_DST_G4HIT_D3mV_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.hitmap", &DST_G4HIT_D3mV_hitmap, &b_DST_G4HIT_D3mV_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.layers", &DST_G4HIT_D3mV_layers, &b_DST_G4HIT_D3mV_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mV.layerMaxID", &DST_G4HIT_D3mV_layerMaxID, &b_DST_G4HIT_D3mV_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.fUniqueID", &DST_G4HIT_D3mXp_fUniqueID, &b_DST_G4HIT_D3mXp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.fBits", &DST_G4HIT_D3mXp_fBits, &b_DST_G4HIT_D3mXp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.id", &DST_G4HIT_D3mXp_id, &b_DST_G4HIT_D3mXp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.hitmap", &DST_G4HIT_D3mXp_hitmap, &b_DST_G4HIT_D3mXp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.layers", &DST_G4HIT_D3mXp_layers, &b_DST_G4HIT_D3mXp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mXp.layerMaxID", &DST_G4HIT_D3mXp_layerMaxID, &b_DST_G4HIT_D3mXp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.fUniqueID", &DST_G4HIT_D3mX_fUniqueID, &b_DST_G4HIT_D3mX_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.fBits", &DST_G4HIT_D3mX_fBits, &b_DST_G4HIT_D3mX_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.id", &DST_G4HIT_D3mX_id, &b_DST_G4HIT_D3mX_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.hitmap", &DST_G4HIT_D3mX_hitmap, &b_DST_G4HIT_D3mX_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.layers", &DST_G4HIT_D3mX_layers, &b_DST_G4HIT_D3mX_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mX.layerMaxID", &DST_G4HIT_D3mX_layerMaxID, &b_DST_G4HIT_D3mX_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.fUniqueID", &DST_G4HIT_D3mUp_fUniqueID, &b_DST_G4HIT_D3mUp_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.fBits", &DST_G4HIT_D3mUp_fBits, &b_DST_G4HIT_D3mUp_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.id", &DST_G4HIT_D3mUp_id, &b_DST_G4HIT_D3mUp_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.hitmap", &DST_G4HIT_D3mUp_hitmap, &b_DST_G4HIT_D3mUp_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.layers", &DST_G4HIT_D3mUp_layers, &b_DST_G4HIT_D3mUp_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mUp.layerMaxID", &DST_G4HIT_D3mUp_layerMaxID, &b_DST_G4HIT_D3mUp_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.fUniqueID", &DST_G4HIT_D3mU_fUniqueID, &b_DST_G4HIT_D3mU_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.fBits", &DST_G4HIT_D3mU_fBits, &b_DST_G4HIT_D3mU_fBits);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.id", &DST_G4HIT_D3mU_id, &b_DST_G4HIT_D3mU_id);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.hitmap", &DST_G4HIT_D3mU_hitmap, &b_DST_G4HIT_D3mU_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.layers", &DST_G4HIT_D3mU_layers, &b_DST_G4HIT_D3mU_layers);
   fChain->SetBranchAddress("DST.G4HIT_D3mU.layerMaxID", &DST_G4HIT_D3mU_layerMaxID, &b_DST_G4HIT_D3mU_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H1B.fUniqueID", &DST_G4HIT_H1B_fUniqueID, &b_DST_G4HIT_H1B_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H1B.fBits", &DST_G4HIT_H1B_fBits, &b_DST_G4HIT_H1B_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H1B.id", &DST_G4HIT_H1B_id, &b_DST_G4HIT_H1B_id);
   fChain->SetBranchAddress("DST.G4HIT_H1B.hitmap", &DST_G4HIT_H1B_hitmap, &b_DST_G4HIT_H1B_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H1B.layers", &DST_G4HIT_H1B_layers, &b_DST_G4HIT_H1B_layers);
   fChain->SetBranchAddress("DST.G4HIT_H1B.layerMaxID", &DST_G4HIT_H1B_layerMaxID, &b_DST_G4HIT_H1B_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H1T.fUniqueID", &DST_G4HIT_H1T_fUniqueID, &b_DST_G4HIT_H1T_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H1T.fBits", &DST_G4HIT_H1T_fBits, &b_DST_G4HIT_H1T_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H1T.id", &DST_G4HIT_H1T_id, &b_DST_G4HIT_H1T_id);
   fChain->SetBranchAddress("DST.G4HIT_H1T.hitmap", &DST_G4HIT_H1T_hitmap, &b_DST_G4HIT_H1T_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H1T.layers", &DST_G4HIT_H1T_layers, &b_DST_G4HIT_H1T_layers);
   fChain->SetBranchAddress("DST.G4HIT_H1T.layerMaxID", &DST_G4HIT_H1T_layerMaxID, &b_DST_G4HIT_H1T_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H1L.fUniqueID", &DST_G4HIT_H1L_fUniqueID, &b_DST_G4HIT_H1L_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H1L.fBits", &DST_G4HIT_H1L_fBits, &b_DST_G4HIT_H1L_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H1L.id", &DST_G4HIT_H1L_id, &b_DST_G4HIT_H1L_id);
   fChain->SetBranchAddress("DST.G4HIT_H1L.hitmap", &DST_G4HIT_H1L_hitmap, &b_DST_G4HIT_H1L_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H1L.layers", &DST_G4HIT_H1L_layers, &b_DST_G4HIT_H1L_layers);
   fChain->SetBranchAddress("DST.G4HIT_H1L.layerMaxID", &DST_G4HIT_H1L_layerMaxID, &b_DST_G4HIT_H1L_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H1R.fUniqueID", &DST_G4HIT_H1R_fUniqueID, &b_DST_G4HIT_H1R_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H1R.fBits", &DST_G4HIT_H1R_fBits, &b_DST_G4HIT_H1R_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H1R.id", &DST_G4HIT_H1R_id, &b_DST_G4HIT_H1R_id);
   fChain->SetBranchAddress("DST.G4HIT_H1R.hitmap", &DST_G4HIT_H1R_hitmap, &b_DST_G4HIT_H1R_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H1R.layers", &DST_G4HIT_H1R_layers, &b_DST_G4HIT_H1R_layers);
   fChain->SetBranchAddress("DST.G4HIT_H1R.layerMaxID", &DST_G4HIT_H1R_layerMaxID, &b_DST_G4HIT_H1R_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H2L.fUniqueID", &DST_G4HIT_H2L_fUniqueID, &b_DST_G4HIT_H2L_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H2L.fBits", &DST_G4HIT_H2L_fBits, &b_DST_G4HIT_H2L_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H2L.id", &DST_G4HIT_H2L_id, &b_DST_G4HIT_H2L_id);
   fChain->SetBranchAddress("DST.G4HIT_H2L.hitmap", &DST_G4HIT_H2L_hitmap, &b_DST_G4HIT_H2L_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H2L.layers", &DST_G4HIT_H2L_layers, &b_DST_G4HIT_H2L_layers);
   fChain->SetBranchAddress("DST.G4HIT_H2L.layerMaxID", &DST_G4HIT_H2L_layerMaxID, &b_DST_G4HIT_H2L_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H2R.fUniqueID", &DST_G4HIT_H2R_fUniqueID, &b_DST_G4HIT_H2R_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H2R.fBits", &DST_G4HIT_H2R_fBits, &b_DST_G4HIT_H2R_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H2R.id", &DST_G4HIT_H2R_id, &b_DST_G4HIT_H2R_id);
   fChain->SetBranchAddress("DST.G4HIT_H2R.hitmap", &DST_G4HIT_H2R_hitmap, &b_DST_G4HIT_H2R_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H2R.layers", &DST_G4HIT_H2R_layers, &b_DST_G4HIT_H2R_layers);
   fChain->SetBranchAddress("DST.G4HIT_H2R.layerMaxID", &DST_G4HIT_H2R_layerMaxID, &b_DST_G4HIT_H2R_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H2B.fUniqueID", &DST_G4HIT_H2B_fUniqueID, &b_DST_G4HIT_H2B_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H2B.fBits", &DST_G4HIT_H2B_fBits, &b_DST_G4HIT_H2B_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H2B.id", &DST_G4HIT_H2B_id, &b_DST_G4HIT_H2B_id);
   fChain->SetBranchAddress("DST.G4HIT_H2B.hitmap", &DST_G4HIT_H2B_hitmap, &b_DST_G4HIT_H2B_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H2B.layers", &DST_G4HIT_H2B_layers, &b_DST_G4HIT_H2B_layers);
   fChain->SetBranchAddress("DST.G4HIT_H2B.layerMaxID", &DST_G4HIT_H2B_layerMaxID, &b_DST_G4HIT_H2B_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H2T.fUniqueID", &DST_G4HIT_H2T_fUniqueID, &b_DST_G4HIT_H2T_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H2T.fBits", &DST_G4HIT_H2T_fBits, &b_DST_G4HIT_H2T_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H2T.id", &DST_G4HIT_H2T_id, &b_DST_G4HIT_H2T_id);
   fChain->SetBranchAddress("DST.G4HIT_H2T.hitmap", &DST_G4HIT_H2T_hitmap, &b_DST_G4HIT_H2T_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H2T.layers", &DST_G4HIT_H2T_layers, &b_DST_G4HIT_H2T_layers);
   fChain->SetBranchAddress("DST.G4HIT_H2T.layerMaxID", &DST_G4HIT_H2T_layerMaxID, &b_DST_G4HIT_H2T_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H3B.fUniqueID", &DST_G4HIT_H3B_fUniqueID, &b_DST_G4HIT_H3B_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H3B.fBits", &DST_G4HIT_H3B_fBits, &b_DST_G4HIT_H3B_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H3B.id", &DST_G4HIT_H3B_id, &b_DST_G4HIT_H3B_id);
   fChain->SetBranchAddress("DST.G4HIT_H3B.hitmap", &DST_G4HIT_H3B_hitmap, &b_DST_G4HIT_H3B_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H3B.layers", &DST_G4HIT_H3B_layers, &b_DST_G4HIT_H3B_layers);
   fChain->SetBranchAddress("DST.G4HIT_H3B.layerMaxID", &DST_G4HIT_H3B_layerMaxID, &b_DST_G4HIT_H3B_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H3T.fUniqueID", &DST_G4HIT_H3T_fUniqueID, &b_DST_G4HIT_H3T_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H3T.fBits", &DST_G4HIT_H3T_fBits, &b_DST_G4HIT_H3T_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H3T.id", &DST_G4HIT_H3T_id, &b_DST_G4HIT_H3T_id);
   fChain->SetBranchAddress("DST.G4HIT_H3T.hitmap", &DST_G4HIT_H3T_hitmap, &b_DST_G4HIT_H3T_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H3T.layers", &DST_G4HIT_H3T_layers, &b_DST_G4HIT_H3T_layers);
   fChain->SetBranchAddress("DST.G4HIT_H3T.layerMaxID", &DST_G4HIT_H3T_layerMaxID, &b_DST_G4HIT_H3T_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.fUniqueID", &DST_G4HIT_H4Y1L_fUniqueID, &b_DST_G4HIT_H4Y1L_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.fBits", &DST_G4HIT_H4Y1L_fBits, &b_DST_G4HIT_H4Y1L_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.id", &DST_G4HIT_H4Y1L_id, &b_DST_G4HIT_H4Y1L_id);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.hitmap", &DST_G4HIT_H4Y1L_hitmap, &b_DST_G4HIT_H4Y1L_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.layers", &DST_G4HIT_H4Y1L_layers, &b_DST_G4HIT_H4Y1L_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1L.layerMaxID", &DST_G4HIT_H4Y1L_layerMaxID, &b_DST_G4HIT_H4Y1L_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.fUniqueID", &DST_G4HIT_H4Y1R_fUniqueID, &b_DST_G4HIT_H4Y1R_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.fBits", &DST_G4HIT_H4Y1R_fBits, &b_DST_G4HIT_H4Y1R_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.id", &DST_G4HIT_H4Y1R_id, &b_DST_G4HIT_H4Y1R_id);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.hitmap", &DST_G4HIT_H4Y1R_hitmap, &b_DST_G4HIT_H4Y1R_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.layers", &DST_G4HIT_H4Y1R_layers, &b_DST_G4HIT_H4Y1R_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4Y1R.layerMaxID", &DST_G4HIT_H4Y1R_layerMaxID, &b_DST_G4HIT_H4Y1R_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.fUniqueID", &DST_G4HIT_H4Y2L_fUniqueID, &b_DST_G4HIT_H4Y2L_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.fBits", &DST_G4HIT_H4Y2L_fBits, &b_DST_G4HIT_H4Y2L_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.id", &DST_G4HIT_H4Y2L_id, &b_DST_G4HIT_H4Y2L_id);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.hitmap", &DST_G4HIT_H4Y2L_hitmap, &b_DST_G4HIT_H4Y2L_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.layers", &DST_G4HIT_H4Y2L_layers, &b_DST_G4HIT_H4Y2L_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2L.layerMaxID", &DST_G4HIT_H4Y2L_layerMaxID, &b_DST_G4HIT_H4Y2L_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.fUniqueID", &DST_G4HIT_H4Y2R_fUniqueID, &b_DST_G4HIT_H4Y2R_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.fBits", &DST_G4HIT_H4Y2R_fBits, &b_DST_G4HIT_H4Y2R_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.id", &DST_G4HIT_H4Y2R_id, &b_DST_G4HIT_H4Y2R_id);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.hitmap", &DST_G4HIT_H4Y2R_hitmap, &b_DST_G4HIT_H4Y2R_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.layers", &DST_G4HIT_H4Y2R_layers, &b_DST_G4HIT_H4Y2R_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4Y2R.layerMaxID", &DST_G4HIT_H4Y2R_layerMaxID, &b_DST_G4HIT_H4Y2R_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4B.fUniqueID", &DST_G4HIT_H4B_fUniqueID, &b_DST_G4HIT_H4B_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4B.fBits", &DST_G4HIT_H4B_fBits, &b_DST_G4HIT_H4B_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4B.id", &DST_G4HIT_H4B_id, &b_DST_G4HIT_H4B_id);
   fChain->SetBranchAddress("DST.G4HIT_H4B.hitmap", &DST_G4HIT_H4B_hitmap, &b_DST_G4HIT_H4B_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4B.layers", &DST_G4HIT_H4B_layers, &b_DST_G4HIT_H4B_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4B.layerMaxID", &DST_G4HIT_H4B_layerMaxID, &b_DST_G4HIT_H4B_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_H4T.fUniqueID", &DST_G4HIT_H4T_fUniqueID, &b_DST_G4HIT_H4T_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_H4T.fBits", &DST_G4HIT_H4T_fBits, &b_DST_G4HIT_H4T_fBits);
   fChain->SetBranchAddress("DST.G4HIT_H4T.id", &DST_G4HIT_H4T_id, &b_DST_G4HIT_H4T_id);
   fChain->SetBranchAddress("DST.G4HIT_H4T.hitmap", &DST_G4HIT_H4T_hitmap, &b_DST_G4HIT_H4T_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_H4T.layers", &DST_G4HIT_H4T_layers, &b_DST_G4HIT_H4T_layers);
   fChain->SetBranchAddress("DST.G4HIT_H4T.layerMaxID", &DST_G4HIT_H4T_layerMaxID, &b_DST_G4HIT_H4T_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.fUniqueID", &DST_G4HIT_P1Y1_fUniqueID, &b_DST_G4HIT_P1Y1_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.fBits", &DST_G4HIT_P1Y1_fBits, &b_DST_G4HIT_P1Y1_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.id", &DST_G4HIT_P1Y1_id, &b_DST_G4HIT_P1Y1_id);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.hitmap", &DST_G4HIT_P1Y1_hitmap, &b_DST_G4HIT_P1Y1_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.layers", &DST_G4HIT_P1Y1_layers, &b_DST_G4HIT_P1Y1_layers);
   fChain->SetBranchAddress("DST.G4HIT_P1Y1.layerMaxID", &DST_G4HIT_P1Y1_layerMaxID, &b_DST_G4HIT_P1Y1_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.fUniqueID", &DST_G4HIT_P1Y2_fUniqueID, &b_DST_G4HIT_P1Y2_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.fBits", &DST_G4HIT_P1Y2_fBits, &b_DST_G4HIT_P1Y2_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.id", &DST_G4HIT_P1Y2_id, &b_DST_G4HIT_P1Y2_id);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.hitmap", &DST_G4HIT_P1Y2_hitmap, &b_DST_G4HIT_P1Y2_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.layers", &DST_G4HIT_P1Y2_layers, &b_DST_G4HIT_P1Y2_layers);
   fChain->SetBranchAddress("DST.G4HIT_P1Y2.layerMaxID", &DST_G4HIT_P1Y2_layerMaxID, &b_DST_G4HIT_P1Y2_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.fUniqueID", &DST_G4HIT_P1X1_fUniqueID, &b_DST_G4HIT_P1X1_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.fBits", &DST_G4HIT_P1X1_fBits, &b_DST_G4HIT_P1X1_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.id", &DST_G4HIT_P1X1_id, &b_DST_G4HIT_P1X1_id);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.hitmap", &DST_G4HIT_P1X1_hitmap, &b_DST_G4HIT_P1X1_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.layers", &DST_G4HIT_P1X1_layers, &b_DST_G4HIT_P1X1_layers);
   fChain->SetBranchAddress("DST.G4HIT_P1X1.layerMaxID", &DST_G4HIT_P1X1_layerMaxID, &b_DST_G4HIT_P1X1_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.fUniqueID", &DST_G4HIT_P1X2_fUniqueID, &b_DST_G4HIT_P1X2_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.fBits", &DST_G4HIT_P1X2_fBits, &b_DST_G4HIT_P1X2_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.id", &DST_G4HIT_P1X2_id, &b_DST_G4HIT_P1X2_id);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.hitmap", &DST_G4HIT_P1X2_hitmap, &b_DST_G4HIT_P1X2_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.layers", &DST_G4HIT_P1X2_layers, &b_DST_G4HIT_P1X2_layers);
   fChain->SetBranchAddress("DST.G4HIT_P1X2.layerMaxID", &DST_G4HIT_P1X2_layerMaxID, &b_DST_G4HIT_P1X2_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.fUniqueID", &DST_G4HIT_P2X1_fUniqueID, &b_DST_G4HIT_P2X1_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.fBits", &DST_G4HIT_P2X1_fBits, &b_DST_G4HIT_P2X1_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.id", &DST_G4HIT_P2X1_id, &b_DST_G4HIT_P2X1_id);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.hitmap", &DST_G4HIT_P2X1_hitmap, &b_DST_G4HIT_P2X1_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.layers", &DST_G4HIT_P2X1_layers, &b_DST_G4HIT_P2X1_layers);
   fChain->SetBranchAddress("DST.G4HIT_P2X1.layerMaxID", &DST_G4HIT_P2X1_layerMaxID, &b_DST_G4HIT_P2X1_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.fUniqueID", &DST_G4HIT_P2X2_fUniqueID, &b_DST_G4HIT_P2X2_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.fBits", &DST_G4HIT_P2X2_fBits, &b_DST_G4HIT_P2X2_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.id", &DST_G4HIT_P2X2_id, &b_DST_G4HIT_P2X2_id);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.hitmap", &DST_G4HIT_P2X2_hitmap, &b_DST_G4HIT_P2X2_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.layers", &DST_G4HIT_P2X2_layers, &b_DST_G4HIT_P2X2_layers);
   fChain->SetBranchAddress("DST.G4HIT_P2X2.layerMaxID", &DST_G4HIT_P2X2_layerMaxID, &b_DST_G4HIT_P2X2_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.fUniqueID", &DST_G4HIT_P2Y1_fUniqueID, &b_DST_G4HIT_P2Y1_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.fBits", &DST_G4HIT_P2Y1_fBits, &b_DST_G4HIT_P2Y1_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.id", &DST_G4HIT_P2Y1_id, &b_DST_G4HIT_P2Y1_id);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.hitmap", &DST_G4HIT_P2Y1_hitmap, &b_DST_G4HIT_P2Y1_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.layers", &DST_G4HIT_P2Y1_layers, &b_DST_G4HIT_P2Y1_layers);
   fChain->SetBranchAddress("DST.G4HIT_P2Y1.layerMaxID", &DST_G4HIT_P2Y1_layerMaxID, &b_DST_G4HIT_P2Y1_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.fUniqueID", &DST_G4HIT_P2Y2_fUniqueID, &b_DST_G4HIT_P2Y2_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.fBits", &DST_G4HIT_P2Y2_fBits, &b_DST_G4HIT_P2Y2_fBits);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.id", &DST_G4HIT_P2Y2_id, &b_DST_G4HIT_P2Y2_id);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.hitmap", &DST_G4HIT_P2Y2_hitmap, &b_DST_G4HIT_P2Y2_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.layers", &DST_G4HIT_P2Y2_layers, &b_DST_G4HIT_P2Y2_layers);
   fChain->SetBranchAddress("DST.G4HIT_P2Y2.layerMaxID", &DST_G4HIT_P2Y2_layerMaxID, &b_DST_G4HIT_P2Y2_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.fUniqueID", &DST_G4HIT_DP1TL_fUniqueID, &b_DST_G4HIT_DP1TL_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.fBits", &DST_G4HIT_DP1TL_fBits, &b_DST_G4HIT_DP1TL_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.id", &DST_G4HIT_DP1TL_id, &b_DST_G4HIT_DP1TL_id);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.hitmap", &DST_G4HIT_DP1TL_hitmap, &b_DST_G4HIT_DP1TL_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.layers", &DST_G4HIT_DP1TL_layers, &b_DST_G4HIT_DP1TL_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP1TL.layerMaxID", &DST_G4HIT_DP1TL_layerMaxID, &b_DST_G4HIT_DP1TL_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.fUniqueID", &DST_G4HIT_DP1TR_fUniqueID, &b_DST_G4HIT_DP1TR_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.fBits", &DST_G4HIT_DP1TR_fBits, &b_DST_G4HIT_DP1TR_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.id", &DST_G4HIT_DP1TR_id, &b_DST_G4HIT_DP1TR_id);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.hitmap", &DST_G4HIT_DP1TR_hitmap, &b_DST_G4HIT_DP1TR_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.layers", &DST_G4HIT_DP1TR_layers, &b_DST_G4HIT_DP1TR_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP1TR.layerMaxID", &DST_G4HIT_DP1TR_layerMaxID, &b_DST_G4HIT_DP1TR_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.fUniqueID", &DST_G4HIT_DP1BL_fUniqueID, &b_DST_G4HIT_DP1BL_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.fBits", &DST_G4HIT_DP1BL_fBits, &b_DST_G4HIT_DP1BL_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.id", &DST_G4HIT_DP1BL_id, &b_DST_G4HIT_DP1BL_id);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.hitmap", &DST_G4HIT_DP1BL_hitmap, &b_DST_G4HIT_DP1BL_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.layers", &DST_G4HIT_DP1BL_layers, &b_DST_G4HIT_DP1BL_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP1BL.layerMaxID", &DST_G4HIT_DP1BL_layerMaxID, &b_DST_G4HIT_DP1BL_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.fUniqueID", &DST_G4HIT_DP1BR_fUniqueID, &b_DST_G4HIT_DP1BR_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.fBits", &DST_G4HIT_DP1BR_fBits, &b_DST_G4HIT_DP1BR_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.id", &DST_G4HIT_DP1BR_id, &b_DST_G4HIT_DP1BR_id);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.hitmap", &DST_G4HIT_DP1BR_hitmap, &b_DST_G4HIT_DP1BR_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.layers", &DST_G4HIT_DP1BR_layers, &b_DST_G4HIT_DP1BR_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP1BR.layerMaxID", &DST_G4HIT_DP1BR_layerMaxID, &b_DST_G4HIT_DP1BR_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.fUniqueID", &DST_G4HIT_DP2TL_fUniqueID, &b_DST_G4HIT_DP2TL_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.fBits", &DST_G4HIT_DP2TL_fBits, &b_DST_G4HIT_DP2TL_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.id", &DST_G4HIT_DP2TL_id, &b_DST_G4HIT_DP2TL_id);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.hitmap", &DST_G4HIT_DP2TL_hitmap, &b_DST_G4HIT_DP2TL_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.layers", &DST_G4HIT_DP2TL_layers, &b_DST_G4HIT_DP2TL_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP2TL.layerMaxID", &DST_G4HIT_DP2TL_layerMaxID, &b_DST_G4HIT_DP2TL_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.fUniqueID", &DST_G4HIT_DP2TR_fUniqueID, &b_DST_G4HIT_DP2TR_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.fBits", &DST_G4HIT_DP2TR_fBits, &b_DST_G4HIT_DP2TR_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.id", &DST_G4HIT_DP2TR_id, &b_DST_G4HIT_DP2TR_id);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.hitmap", &DST_G4HIT_DP2TR_hitmap, &b_DST_G4HIT_DP2TR_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.layers", &DST_G4HIT_DP2TR_layers, &b_DST_G4HIT_DP2TR_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP2TR.layerMaxID", &DST_G4HIT_DP2TR_layerMaxID, &b_DST_G4HIT_DP2TR_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.fUniqueID", &DST_G4HIT_DP2BL_fUniqueID, &b_DST_G4HIT_DP2BL_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.fBits", &DST_G4HIT_DP2BL_fBits, &b_DST_G4HIT_DP2BL_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.id", &DST_G4HIT_DP2BL_id, &b_DST_G4HIT_DP2BL_id);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.hitmap", &DST_G4HIT_DP2BL_hitmap, &b_DST_G4HIT_DP2BL_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.layers", &DST_G4HIT_DP2BL_layers, &b_DST_G4HIT_DP2BL_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP2BL.layerMaxID", &DST_G4HIT_DP2BL_layerMaxID, &b_DST_G4HIT_DP2BL_layerMaxID);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.fUniqueID", &DST_G4HIT_DP2BR_fUniqueID, &b_DST_G4HIT_DP2BR_fUniqueID);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.fBits", &DST_G4HIT_DP2BR_fBits, &b_DST_G4HIT_DP2BR_fBits);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.id", &DST_G4HIT_DP2BR_id, &b_DST_G4HIT_DP2BR_id);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.hitmap", &DST_G4HIT_DP2BR_hitmap, &b_DST_G4HIT_DP2BR_hitmap);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.layers", &DST_G4HIT_DP2BR_layers, &b_DST_G4HIT_DP2BR_layers);
   fChain->SetBranchAddress("DST.G4HIT_DP2BR.layerMaxID", &DST_G4HIT_DP2BR_layerMaxID, &b_DST_G4HIT_DP2BR_layerMaxID);
   fChain->SetBranchAddress("DST.G4TruthInfo.fUniqueID", &DST_G4TruthInfo_fUniqueID, &b_DST_G4TruthInfo_fUniqueID);
   fChain->SetBranchAddress("DST.G4TruthInfo.fBits", &DST_G4TruthInfo_fBits, &b_DST_G4TruthInfo_fBits);
   fChain->SetBranchAddress("DST.G4TruthInfo.particlemap", &DST_G4TruthInfo_particlemap, &b_DST_G4TruthInfo_particlemap);
   fChain->SetBranchAddress("DST.G4TruthInfo.vtxmap", &DST_G4TruthInfo_vtxmap, &b_DST_G4TruthInfo_vtxmap);
   fChain->SetBranchAddress("DST.G4TruthInfo.showermap", &DST_G4TruthInfo_showermap, &b_DST_G4TruthInfo_showermap);
   fChain->SetBranchAddress("DST.G4TruthInfo.particle_embed_flags", &DST_G4TruthInfo_particle_embed_flags, &b_DST_G4TruthInfo_particle_embed_flags);
   fChain->SetBranchAddress("DST.G4TruthInfo.vertex_embed_flags", &DST_G4TruthInfo_vertex_embed_flags, &b_DST_G4TruthInfo_vertex_embed_flags);
   fChain->SetBranchAddress("DST.SQHitVector.fUniqueID", &DST_SQHitVector_fUniqueID, &b_DST_SQHitVector_fUniqueID);
   fChain->SetBranchAddress("DST.SQHitVector.fBits", &DST_SQHitVector_fBits, &b_DST_SQHitVector_fBits);
   fChain->SetBranchAddress("DST.SQHitVector._vector", &DST_SQHitVector__vector, &b_DST_SQHitVector__vector);
   fChain->SetBranchAddress("DST.SQEvent.fUniqueID", &DST_SQEvent_fUniqueID, &b_DST_SQEvent_fUniqueID);
   fChain->SetBranchAddress("DST.SQEvent.fBits", &DST_SQEvent_fBits, &b_DST_SQEvent_fBits);
   fChain->SetBranchAddress("DST.SQEvent._run_id", &DST_SQEvent__run_id, &b_DST_SQEvent__run_id);
   fChain->SetBranchAddress("DST.SQEvent._spill_id", &DST_SQEvent__spill_id, &b_DST_SQEvent__spill_id);
   fChain->SetBranchAddress("DST.SQEvent._event_id", &DST_SQEvent__event_id, &b_DST_SQEvent__event_id);
   fChain->SetBranchAddress("DST.SQEvent._coda_event_id", &DST_SQEvent__coda_event_id, &b_DST_SQEvent__coda_event_id);
   fChain->SetBranchAddress("DST.SQEvent._trigger", &DST_SQEvent__trigger, &b_DST_SQEvent__trigger);
   fChain->SetBranchAddress("DST.SQEvent._raw_matrix[5]", DST_SQEvent__raw_matrix, &b_DST_SQEvent__raw_matrix);
   fChain->SetBranchAddress("DST.SQEvent._after_inh_matrix[5]", DST_SQEvent__after_inh_matrix, &b_DST_SQEvent__after_inh_matrix);
   fChain->SetBranchAddress("DST.SQEvent._data_quality", &DST_SQEvent__data_quality, &b_DST_SQEvent__data_quality);
   fChain->SetBranchAddress("DST.SQEvent._vme_time", &DST_SQEvent__vme_time, &b_DST_SQEvent__vme_time);
   fChain->SetBranchAddress("DST.SQEvent._qie_presums[4]", DST_SQEvent__qie_presums, &b_DST_SQEvent__qie_presums);
   fChain->SetBranchAddress("DST.SQEvent._qie_trig_cnt", &DST_SQEvent__qie_trig_cnt, &b_DST_SQEvent__qie_trig_cnt);
   fChain->SetBranchAddress("DST.SQEvent._qie_turn_id", &DST_SQEvent__qie_turn_id, &b_DST_SQEvent__qie_turn_id);
   fChain->SetBranchAddress("DST.SQEvent._qie_rf_id", &DST_SQEvent__qie_rf_id, &b_DST_SQEvent__qie_rf_id);
   fChain->SetBranchAddress("DST.SQEvent._qie_rf_inte[33]", DST_SQEvent__qie_rf_inte, &b_DST_SQEvent__qie_rf_inte);
   fChain->SetBranchAddress("DST.SQEvent._flag_v1495", &DST_SQEvent__flag_v1495, &b_DST_SQEvent__flag_v1495);
   fChain->SetBranchAddress("DST.SQEvent._n_board_qie", &DST_SQEvent__n_board_qie, &b_DST_SQEvent__n_board_qie);
   fChain->SetBranchAddress("DST.SQEvent._n_board_v1495", &DST_SQEvent__n_board_v1495, &b_DST_SQEvent__n_board_v1495);
   fChain->SetBranchAddress("DST.SQEvent._n_board_taiwan", &DST_SQEvent__n_board_taiwan, &b_DST_SQEvent__n_board_taiwan);
   fChain->SetBranchAddress("DST.SQEvent._n_board_trig_b", &DST_SQEvent__n_board_trig_b, &b_DST_SQEvent__n_board_trig_b);
   fChain->SetBranchAddress("DST.SQEvent._n_board_trig_c", &DST_SQEvent__n_board_trig_c, &b_DST_SQEvent__n_board_trig_c);
   fChain->SetBranchAddress("DST.SQMCEvent.fUniqueID", &DST_SQMCEvent_fUniqueID, &b_DST_SQMCEvent_fUniqueID);
   fChain->SetBranchAddress("DST.SQMCEvent.fBits", &DST_SQMCEvent_fBits, &b_DST_SQMCEvent_fBits);
   fChain->SetBranchAddress("DST.SQMCEvent._proc_id", &DST_SQMCEvent__proc_id, &b_DST_SQMCEvent__proc_id);
   fChain->SetBranchAddress("DST.SQMCEvent._xsec", &DST_SQMCEvent__xsec, &b_DST_SQMCEvent__xsec);
   fChain->SetBranchAddress("DST.SQMCEvent._weight", &DST_SQMCEvent__weight, &b_DST_SQMCEvent__weight);
   fChain->SetBranchAddress("DST.SQMCEvent._par_id[4]", DST_SQMCEvent__par_id, &b_DST_SQMCEvent__par_id);
   fChain->SetBranchAddress("DST.SQMCEvent._par_mom[4]", DST_SQMCEvent__par_mom, &b_DST_SQMCEvent__par_mom);
   fChain->SetBranchAddress("DST.SQTruthTrackVector.fUniqueID", &DST_SQTruthTrackVector_fUniqueID, &b_DST_SQTruthTrackVector_fUniqueID);
   fChain->SetBranchAddress("DST.SQTruthTrackVector.fBits", &DST_SQTruthTrackVector_fBits, &b_DST_SQTruthTrackVector_fBits);
   fChain->SetBranchAddress("DST.SQTruthTrackVector._vector", &DST_SQTruthTrackVector__vector, &b_DST_SQTruthTrackVector__vector);
   fChain->SetBranchAddress("DST.SQTruthDimuonVector.fUniqueID", &DST_SQTruthDimuonVector_fUniqueID, &b_DST_SQTruthDimuonVector_fUniqueID);
   fChain->SetBranchAddress("DST.SQTruthDimuonVector.fBits", &DST_SQTruthDimuonVector_fBits, &b_DST_SQTruthDimuonVector_fBits);
   fChain->SetBranchAddress("DST.SQTruthDimuonVector._vector", &DST_SQTruthDimuonVector__vector, &b_DST_SQTruthDimuonVector__vector);
   fChain->SetBranchAddress("DST.SRecEvent.fUniqueID", &DST_SRecEvent_fUniqueID, &b_DST_SRecEvent_fUniqueID);
   fChain->SetBranchAddress("DST.SRecEvent.fBits", &DST_SRecEvent_fBits, &b_DST_SRecEvent_fBits);
   fChain->SetBranchAddress("DST.SRecEvent.fRecStatus", &DST_SRecEvent_fRecStatus, &b_DST_SRecEvent_fRecStatus);
   fChain->SetBranchAddress("DST.SRecEvent.fRunID", &DST_SRecEvent_fRunID, &b_DST_SRecEvent_fRunID);
   fChain->SetBranchAddress("DST.SRecEvent.fSpillID", &DST_SRecEvent_fSpillID, &b_DST_SRecEvent_fSpillID);
   fChain->SetBranchAddress("DST.SRecEvent.fEventID", &DST_SRecEvent_fEventID, &b_DST_SRecEvent_fEventID);
   fChain->SetBranchAddress("DST.SRecEvent.fTargetPos", &DST_SRecEvent_fTargetPos, &b_DST_SRecEvent_fTargetPos);
   fChain->SetBranchAddress("DST.SRecEvent.fTriggerBits", &DST_SRecEvent_fTriggerBits, &b_DST_SRecEvent_fTriggerBits);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks", &DST_SRecEvent_fAllTracks_, &b_DST_SRecEvent_fAllTracks_);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fUniqueID", DST_SRecEvent_fAllTracks_fUniqueID, &b_DST_SRecEvent_fAllTracks_fUniqueID);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fBits", DST_SRecEvent_fAllTracks_fBits, &b_DST_SRecEvent_fAllTracks_fBits);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisq", DST_SRecEvent_fAllTracks_fChisq, &b_DST_SRecEvent_fAllTracks_fChisq);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fHitIndex", DST_SRecEvent_fAllTracks_fHitIndex, &b_DST_SRecEvent_fAllTracks_fHitIndex);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fState", DST_SRecEvent_fAllTracks_fState, &b_DST_SRecEvent_fAllTracks_fState);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fCovar", DST_SRecEvent_fAllTracks_fCovar, &b_DST_SRecEvent_fAllTracks_fCovar);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fZ", DST_SRecEvent_fAllTracks_fZ, &b_DST_SRecEvent_fAllTracks_fZ);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisqAtNode", DST_SRecEvent_fAllTracks_fChisqAtNode, &b_DST_SRecEvent_fAllTracks_fChisqAtNode);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fDumpFacePos", DST_SRecEvent_fAllTracks_fDumpFacePos, &b_DST_SRecEvent_fAllTracks_fDumpFacePos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fDumpPos", DST_SRecEvent_fAllTracks_fDumpPos, &b_DST_SRecEvent_fAllTracks_fDumpPos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fTargetPos", DST_SRecEvent_fAllTracks_fTargetPos, &b_DST_SRecEvent_fAllTracks_fTargetPos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fXVertexPos", DST_SRecEvent_fAllTracks_fXVertexPos, &b_DST_SRecEvent_fAllTracks_fXVertexPos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fYVertexPos", DST_SRecEvent_fAllTracks_fYVertexPos, &b_DST_SRecEvent_fAllTracks_fYVertexPos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fDumpFaceMom", DST_SRecEvent_fAllTracks_fDumpFaceMom, &b_DST_SRecEvent_fAllTracks_fDumpFaceMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fDumpMom", DST_SRecEvent_fAllTracks_fDumpMom, &b_DST_SRecEvent_fAllTracks_fDumpMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fTargetMom", DST_SRecEvent_fAllTracks_fTargetMom, &b_DST_SRecEvent_fAllTracks_fTargetMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fXVertexMom", DST_SRecEvent_fAllTracks_fXVertexMom, &b_DST_SRecEvent_fAllTracks_fXVertexMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fYVertexMom", DST_SRecEvent_fAllTracks_fYVertexMom, &b_DST_SRecEvent_fAllTracks_fYVertexMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fVertexMom", DST_SRecEvent_fAllTracks_fVertexMom, &b_DST_SRecEvent_fAllTracks_fVertexMom);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fVertexPos", DST_SRecEvent_fAllTracks_fVertexPos, &b_DST_SRecEvent_fAllTracks_fVertexPos);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisqVertex", DST_SRecEvent_fAllTracks_fChisqVertex, &b_DST_SRecEvent_fAllTracks_fChisqVertex);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fStateVertex", DST_SRecEvent_fAllTracks_fStateVertex, &b_DST_SRecEvent_fAllTracks_fStateVertex);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fCovarVertex", DST_SRecEvent_fAllTracks_fCovarVertex, &b_DST_SRecEvent_fAllTracks_fCovarVertex);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fKalmanStatus", DST_SRecEvent_fAllTracks_fKalmanStatus, &b_DST_SRecEvent_fAllTracks_fKalmanStatus);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fTriggerID", DST_SRecEvent_fAllTracks_fTriggerID, &b_DST_SRecEvent_fAllTracks_fTriggerID);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fNPropHitsX", DST_SRecEvent_fAllTracks_fNPropHitsX, &b_DST_SRecEvent_fAllTracks_fNPropHitsX);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fNPropHitsY", DST_SRecEvent_fAllTracks_fNPropHitsY, &b_DST_SRecEvent_fAllTracks_fNPropHitsY);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fPropSlopeX", DST_SRecEvent_fAllTracks_fPropSlopeX, &b_DST_SRecEvent_fAllTracks_fPropSlopeX);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fPropSlopeY", DST_SRecEvent_fAllTracks_fPropSlopeY, &b_DST_SRecEvent_fAllTracks_fPropSlopeY);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisqTarget", DST_SRecEvent_fAllTracks_fChisqTarget, &b_DST_SRecEvent_fAllTracks_fChisqTarget);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisqDump", DST_SRecEvent_fAllTracks_fChisqDump, &b_DST_SRecEvent_fAllTracks_fChisqDump);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fChisqUpstream", DST_SRecEvent_fAllTracks_fChisqUpstream, &b_DST_SRecEvent_fAllTracks_fChisqUpstream);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fGFDetPlaneVec[3]", DST_SRecEvent_fAllTracks_fGFDetPlaneVec, &b_DST_SRecEvent_fAllTracks_fGFDetPlaneVec);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fGFAuxInfo", DST_SRecEvent_fAllTracks_fGFAuxInfo, &b_DST_SRecEvent_fAllTracks_fGFAuxInfo);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fGFStateVec", DST_SRecEvent_fAllTracks_fGFStateVec, &b_DST_SRecEvent_fAllTracks_fGFStateVec);
   fChain->SetBranchAddress("DST.SRecEvent.fAllTracks.fGFCov", DST_SRecEvent_fAllTracks_fGFCov, &b_DST_SRecEvent_fAllTracks_fGFCov);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons", &DST_SRecEvent_fDimuons_, &b_DST_SRecEvent_fDimuons_);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.fUniqueID", DST_SRecEvent_fDimuons_fUniqueID, &b_DST_SRecEvent_fDimuons_fUniqueID);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.fBits", DST_SRecEvent_fDimuons_fBits, &b_DST_SRecEvent_fDimuons_fBits);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.trackID_pos", DST_SRecEvent_fDimuons_trackID_pos, &b_DST_SRecEvent_fDimuons_trackID_pos);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.trackID_neg", DST_SRecEvent_fDimuons_trackID_neg, &b_DST_SRecEvent_fDimuons_trackID_neg);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.p_pos", DST_SRecEvent_fDimuons_p_pos, &b_DST_SRecEvent_fDimuons_p_pos);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.p_neg", DST_SRecEvent_fDimuons_p_neg, &b_DST_SRecEvent_fDimuons_p_neg);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.p_pos_single", DST_SRecEvent_fDimuons_p_pos_single, &b_DST_SRecEvent_fDimuons_p_pos_single);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.p_neg_single", DST_SRecEvent_fDimuons_p_neg_single, &b_DST_SRecEvent_fDimuons_p_neg_single);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.vtx", DST_SRecEvent_fDimuons_vtx, &b_DST_SRecEvent_fDimuons_vtx);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.vtx_pos", DST_SRecEvent_fDimuons_vtx_pos, &b_DST_SRecEvent_fDimuons_vtx_pos);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.vtx_neg", DST_SRecEvent_fDimuons_vtx_neg, &b_DST_SRecEvent_fDimuons_vtx_neg);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.proj_target_pos", DST_SRecEvent_fDimuons_proj_target_pos, &b_DST_SRecEvent_fDimuons_proj_target_pos);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.proj_dump_pos", DST_SRecEvent_fDimuons_proj_dump_pos, &b_DST_SRecEvent_fDimuons_proj_dump_pos);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.proj_target_neg", DST_SRecEvent_fDimuons_proj_target_neg, &b_DST_SRecEvent_fDimuons_proj_target_neg);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.proj_dump_neg", DST_SRecEvent_fDimuons_proj_dump_neg, &b_DST_SRecEvent_fDimuons_proj_dump_neg);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.mass", DST_SRecEvent_fDimuons_mass, &b_DST_SRecEvent_fDimuons_mass);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.pT", DST_SRecEvent_fDimuons_pT, &b_DST_SRecEvent_fDimuons_pT);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.xF", DST_SRecEvent_fDimuons_xF, &b_DST_SRecEvent_fDimuons_xF);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.x1", DST_SRecEvent_fDimuons_x1, &b_DST_SRecEvent_fDimuons_x1);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.x2", DST_SRecEvent_fDimuons_x2, &b_DST_SRecEvent_fDimuons_x2);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.costh", DST_SRecEvent_fDimuons_costh, &b_DST_SRecEvent_fDimuons_costh);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.phi", DST_SRecEvent_fDimuons_phi, &b_DST_SRecEvent_fDimuons_phi);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.mass_single", DST_SRecEvent_fDimuons_mass_single, &b_DST_SRecEvent_fDimuons_mass_single);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_single", DST_SRecEvent_fDimuons_chisq_single, &b_DST_SRecEvent_fDimuons_chisq_single);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_kf", DST_SRecEvent_fDimuons_chisq_kf, &b_DST_SRecEvent_fDimuons_chisq_kf);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_vx", DST_SRecEvent_fDimuons_chisq_vx, &b_DST_SRecEvent_fDimuons_chisq_vx);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_target", DST_SRecEvent_fDimuons_chisq_target, &b_DST_SRecEvent_fDimuons_chisq_target);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_dump", DST_SRecEvent_fDimuons_chisq_dump, &b_DST_SRecEvent_fDimuons_chisq_dump);
   fChain->SetBranchAddress("DST.SRecEvent.fDimuons.chisq_upstream", DST_SRecEvent_fDimuons_chisq_upstream, &b_DST_SRecEvent_fDimuons_chisq_upstream);
   fChain->SetBranchAddress("DST.SRecEvent.fLocalID", &DST_SRecEvent_fLocalID, &b_DST_SRecEvent_fLocalID);
   fChain->SetBranchAddress("DST.SRecEvent.fSource1", &DST_SRecEvent_fSource1, &b_DST_SRecEvent_fSource1);
   fChain->SetBranchAddress("DST.SRecEvent.fSource2", &DST_SRecEvent_fSource2, &b_DST_SRecEvent_fSource2);
   Notify();
}

Bool_t e1039data_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void e1039data_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t e1039data_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef e1039data_tree_cxx
