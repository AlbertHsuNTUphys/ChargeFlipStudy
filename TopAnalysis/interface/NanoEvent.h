#ifndef _nanoevent_h_
#define _nanoevent_h_

#include "TTree.h"

struct NanoEvent_t
{
  NanoEvent_t()
  { 
    n_tight_muon=0; n_loose_muon=0;
    n_tight_ele=0;  n_loose_ele=0;
    n_tight_jet=0; n_bjet_DeepB=0; n_cjet_DeepB_medium=0;
    PrefireWeight=1.; puWeight=1.; genWeight=1.;
  }

  //HADDER
  UInt_t run,luminosityBlock;
  ULong64_t event;

  //FILTER --- FLAG
  Bool_t Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter, Flag_CSCTightHaloFilter, Flag_CSCTightHaloTrkMuUnvetoFilter, Flag_CSCTightHalo2015Filter, Flag_globalTightHalo2016Filter, Flag_globalSuperTightHalo2016Filter, Flag_HcalStripHaloFilter, Flag_hcalLaserEventFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_goodVertices, Flag_eeBadScFilter, Flag_BadPFMuonFilter, Flag_ecalBadCalibFilter;

  //TRIGGER --- HLT
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, HLT_PFMET120_PFMHT120_IDTight, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, HLT_PFHT500_PFMET100_PFMHT100_IDTight, HLT_PFHT700_PFMET85_PFMHT85_IDTight, HLT_PFHT800_PFMET75_PFMHT75_IDTight, HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, HLT_IsoMu27, HLT_Ele32_WPTight_Gsf_L1DoubleEG, HLT_Ele35_WPTight_Gsf, HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t HLT_passEle32WPTight;

  //WEIGHT & SCALE FACTOR
  Float_t genWeight,puWeight,PrefireWeight;
  Float_t Electron_CutBased_LooseID_SF[500],  Electron_CutBased_MediumID_SF[500], Electron_CutBased_TightID_SF[500], Electron_RECO_SF[500], Muon_TightISO_SF[500], Muon_CutBased_TightID_SF[500];

  //BASIC INFORMATION --- MET
  Float_t MET_phi, MET_pt, MET_sumEt;
  
  //BASIC INFORMATION --- GENPARTON
  UInt_t nGenDressedLepton;
  Float_t GenDressedLepton_pt[500], GenDressedLepton_eta[500], GenDressedLepton_phi[500], GenDressedLepton_mass[500];
  Int_t GenDressedLepton_pdgId[500];

  //BASIC INFORMATION --- COMMON REGION
  Int_t n_tight_muon, n_loose_muon, n_tight_ele, n_loose_ele, n_tight_jet, n_bjet_DeepB, n_cjet_DeepB_medium;
  Float_t HT, j1_pt, j1_eta, j1_phi, j1_mass, j2_pt, j2_eta, j2_phi, j2_mass, j3_pt, j3_eta, j3_phi, j3_mass, j4_pt, j4_eta, j4_phi, j4_mass, DeepB_j1_pt, DeepB_j1_eta, DeepB_j1_phi, DeepB_j1_mass, DeepB_j2_pt, DeepB_j2_eta, DeepB_j2_phi, DeepB_j2_mass, DeepC_medium_j1_pt, DeepC_medium_j1_eta, DeepC_medium_j1_phi, DeepC_medium_j1_mass, DeepC_medium_j2_pt, DeepC_medium_j2_eta, DeepC_medium_j2_phi, DeepC_medium_j2_mass;

  //BASIC INFORMATION --- TTC REGION
  Char_t ttc_nl, ttc_jets;
  Int_t ttc_region, ttc_l1_id, ttc_l2_id, ttc_l1_pdgid, ttc_l2_pdgid;
  Float_t ttc_l1_pt, ttc_l1_eta, ttc_l1_phi, ttc_l1_mass, ttc_l2_pt, ttc_l2_eta, ttc_l2_phi, ttc_l2_mass, ttc_mll, ttc_drll, ttc_dphill, ttc_met, ttc_met_phi, ttc_dr_l1j1, ttc_dr_l1j2, ttc_dr_l1j3, ttc_dr_l1j4, ttc_dr_l2j1, ttc_dr_l2j2, ttc_dr_l2j3, ttc_dr_l2j4, ttc_mllj1, ttc_mllj2, ttc_mllj3, ttc_mllj4;

  //BASIC INFORMATION --- WZ REGION
  Int_t WZ_region, WZ_zl1_id, WZ_zl2_id,WZ_wl_id;
  Float_t WZ_zl1_pt, WZ_zl1_eta, WZ_zl1_phi, WZ_zl1_mass, WZ_zl2_pt, WZ_zl2_eta, WZ_zl2_phi, WZ_zl2_mass, WZ_l3_pt, WZ_l3_eta, WZ_l3_phi, WZ_l3_mass, WZ_Z_mass, WZ_Z_pt, WZ_Z_eta, WZ_Z_phi, WZ_met;

  //BASIC INFORMATION --- DY REGION
  Int_t DY_region, DY_l1_id, DY_l2_id, DY_l1_pdgid, DY_l2_pdgid;
  Float_t DY_l1_pt, DY_l1_pt_raw, DY_l1_eta, DY_l1_phi, DY_l1_mass, DY_l2_pt, DY_l2_pt_raw, DY_l2_eta, DY_l2_phi, DY_l2_mass, DY_z_mass, DY_z_mass_raw, DY_z_pt, DY_z_pt_raw, DY_z_eta, DY_z_eta_raw, DY_z_phi, DY_z_phi_raw, DY_drll;
};
void attachToNanoEventTree(TTree *t, NanoEvent_t &ev, bool full, bool is_mc);
#endif
