#include "TopLJets2015/TopAnalysis/interface/NanoEvent.h"

void attachToNanoEventTree(TTree *t, NanoEvent_t &ev, bool full, bool is_mc)
{
  //EVENT HADDER
  t->SetBranchAddress("run",             &ev.run);
  t->SetBranchAddress("event",           &ev.event);
  t->SetBranchAddress("luminosityBlock", &ev.luminosityBlock);

  //FILTER -- FLAG
  t->SetBranchAddress("Flag_HBHENoiseFilter",                    &ev.Flag_HBHENoiseFilter);
  t->SetBranchAddress("Flag_HBHENoiseIsoFilter",                 &ev.Flag_HBHENoiseIsoFilter);
  t->SetBranchAddress("Flag_CSCTightHaloFilter",                 &ev.Flag_CSCTightHaloFilter);
  t->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter",      &ev.Flag_CSCTightHaloTrkMuUnvetoFilter);
  t->SetBranchAddress("Flag_CSCTightHalo2015Filter",             &ev.Flag_CSCTightHalo2015Filter);
  t->SetBranchAddress("Flag_globalTightHalo2016Filter",          &ev.Flag_globalTightHalo2016Filter);
  t->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",     &ev.Flag_globalSuperTightHalo2016Filter);
  t->SetBranchAddress("Flag_HcalStripHaloFilter",                &ev.Flag_HcalStripHaloFilter);
  t->SetBranchAddress("Flag_hcalLaserEventFilter",               &ev.Flag_hcalLaserEventFilter);
  t->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &ev.Flag_EcalDeadCellTriggerPrimitiveFilter);
  t->SetBranchAddress("Flag_goodVertices",                       &ev.Flag_goodVertices );
  t->SetBranchAddress("Flag_eeBadScFilter",                      &ev.Flag_eeBadScFilter);
  t->SetBranchAddress("Flag_BadPFMuonFilter",                    &ev.Flag_BadPFMuonFilter);
  t->SetBranchAddress("Flag_ecalBadCalibFilter",                 &ev.Flag_ecalBadCalibFilter);

  //TRIGGER -- HLT
  t->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          &ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  t->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",             &ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
  t->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight",                      &ev.HLT_PFMET120_PFMHT120_IDTight);
  t->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",              &ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  t->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight",              &ev.HLT_PFHT500_PFMET100_PFMHT100_IDTight);
  t->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight",                &ev.HLT_PFHT700_PFMET85_PFMHT85_IDTight);
  t->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight",                &ev.HLT_PFHT800_PFMET75_PFMHT75_IDTight);
  t->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",          &ev.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
  t->SetBranchAddress("HLT_IsoMu27",                                        &ev.HLT_IsoMu27);
  t->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG",                   &ev.HLT_Ele32_WPTight_Gsf_L1DoubleEG);
  t->SetBranchAddress("HLT_Ele35_WPTight_Gsf",                              &ev.HLT_Ele35_WPTight_Gsf);
  t->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",  &ev.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
  t->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &ev.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
  t->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &ev.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  t->SetBranchAddress("HLT_passEle32WPTight",                               &ev.HLT_passEle32WPTight);

  //WEIGHT & SCALE FACTOR
  if(is_mc){
      t->SetBranchAddress("Electron_CutBased_LooseID_SF",        ev.Electron_CutBased_LooseID_SF);
      t->SetBranchAddress("Electron_CutBased_MediumID_SF",       ev.Electron_CutBased_MediumID_SF);
      t->SetBranchAddress("Electron_CutBased_TightID_SF",        ev.Electron_CutBased_TightID_SF);
      t->SetBranchAddress("Electron_RECO_SF",                    ev.Electron_RECO_SF);
      t->SetBranchAddress("Muon_TightRelIso_TightIDandIPCut_SF", ev.Muon_TightISO_SF);
      t->SetBranchAddress("Muon_CutBased_TightID_SF",            ev.Muon_CutBased_TightID_SF);
      t->SetBranchAddress("PrefireWeight",                       &ev.PrefireWeight);
      t->SetBranchAddress("puWeight",                            &ev.puWeight);
      t->SetBranchAddress("genWeight",                           &ev.genWeight);
  }
  //BASIC INFORMATION --- MET
  t->SetBranchAddress("MET_phi",   &ev.MET_phi);
  t->SetBranchAddress("MET_pt",    &ev.MET_pt);
  t->SetBranchAddress("MET_sumEt", &ev.MET_sumEt);

  //BASIC INFORMATION --- GENPARTON
  if(is_mc){
      t->SetBranchAddress("nGenDressedLepton",      &ev.nGenDressedLepton);
      t->SetBranchAddress("GenDressedLepton_pt",    ev.GenDressedLepton_pt);
      t->SetBranchAddress("GenDressedLepton_eta",   ev.GenDressedLepton_eta);
      t->SetBranchAddress("GenDressedLepton_phi",   ev.GenDressedLepton_phi);
      t->SetBranchAddress("GenDressedLepton_mass",  ev.GenDressedLepton_mass);
      t->SetBranchAddress("GenDressedLepton_pdgId", ev.GenDressedLepton_pdgId);
   }

  //BASIC INFORMATION --- COMMON REGION
  t->SetBranchAddress("n_tight_muon",         &ev.n_tight_muon);
  t->SetBranchAddress("n_loose_muon",         &ev.n_loose_muon);
  t->SetBranchAddress("n_tight_ele",          &ev.n_tight_ele);
  t->SetBranchAddress("n_loose_ele",          &ev.n_loose_ele);
  t->SetBranchAddress("n_tight_jet",          &ev.n_tight_jet);
  t->SetBranchAddress("n_bjet_DeepB",         &ev.n_bjet_DeepB);
  t->SetBranchAddress("n_cjet_DeepB_medium",  &ev.n_cjet_DeepB_medium);
  t->SetBranchAddress("HT",                   &ev.HT);
  t->SetBranchAddress("j1_pt",                &ev.j1_pt);
  t->SetBranchAddress("j1_eta",               &ev.j1_eta);
  t->SetBranchAddress("j1_phi",               &ev.j1_phi);
  t->SetBranchAddress("j1_mass",              &ev.j1_mass);
  t->SetBranchAddress("j2_pt",                &ev.j2_pt);
  t->SetBranchAddress("j2_eta",               &ev.j2_eta);
  t->SetBranchAddress("j2_phi",               &ev.j2_phi);
  t->SetBranchAddress("j2_mass",              &ev.j2_mass);
  t->SetBranchAddress("j3_pt",                &ev.j3_pt);
  t->SetBranchAddress("j3_eta",               &ev.j3_eta);
  t->SetBranchAddress("j3_phi",               &ev.j3_phi);
  t->SetBranchAddress("j3_mass",              &ev.j3_mass);
  t->SetBranchAddress("j4_pt",                &ev.j4_pt);
  t->SetBranchAddress("j4_eta",               &ev.j4_eta);
  t->SetBranchAddress("j4_phi",               &ev.j4_phi);
  t->SetBranchAddress("j4_mass",              &ev.j4_mass);
  t->SetBranchAddress("DeepB_j1_pt",          &ev.DeepB_j1_pt);
  t->SetBranchAddress("DeepB_j1_eta",         &ev.DeepB_j1_eta);
  t->SetBranchAddress("DeepB_j1_phi",         &ev.DeepB_j1_phi);
  t->SetBranchAddress("DeepB_j1_mass",        &ev.DeepB_j1_mass);
  t->SetBranchAddress("DeepB_j2_pt",          &ev.DeepB_j2_pt);
  t->SetBranchAddress("DeepB_j2_eta",         &ev.DeepB_j2_eta);
  t->SetBranchAddress("DeepB_j2_phi",         &ev.DeepB_j2_phi);
  t->SetBranchAddress("DeepB_j2_mass",        &ev.DeepB_j2_mass);
  t->SetBranchAddress("DeepC_medium_j1_pt",   &ev.DeepC_medium_j1_pt);
  t->SetBranchAddress("DeepC_medium_j1_eta",  &ev.DeepC_medium_j1_eta);
  t->SetBranchAddress("DeepC_medium_j1_phi",  &ev.DeepC_medium_j1_phi);
  t->SetBranchAddress("DeepC_medium_j1_mass", &ev.DeepC_medium_j1_mass);
  t->SetBranchAddress("DeepC_medium_j2_pt",   &ev.DeepC_medium_j2_pt);
  t->SetBranchAddress("DeepC_medium_j2_eta",  &ev.DeepC_medium_j2_eta);
  t->SetBranchAddress("DeepC_medium_j2_phi",  &ev.DeepC_medium_j2_phi);
  t->SetBranchAddress("DeepC_medium_j2_mass", &ev.DeepC_medium_j2_mass);

  //BASIC INFORMATION --- TTC REGION
  t->SetBranchAddress("ttc_nl",       &ev.ttc_nl);
  t->SetBranchAddress("ttc_jets",     &ev.ttc_jets);
  t->SetBranchAddress("ttc_region",   &ev.ttc_region);
  t->SetBranchAddress("ttc_l1_id",    &ev.ttc_l1_id);
  t->SetBranchAddress("ttc_l2_id",    &ev.ttc_l2_id);
  t->SetBranchAddress("ttc_l1_pdgid", &ev.ttc_l1_pdgid);
  t->SetBranchAddress("ttc_l2_pdgid", &ev.ttc_l2_pdgid);
  t->SetBranchAddress("ttc_l1_pt",    &ev.ttc_l1_pt);
  t->SetBranchAddress("ttc_l1_eta",   &ev.ttc_l1_eta);
  t->SetBranchAddress("ttc_l1_phi",   &ev.ttc_l1_phi);
  t->SetBranchAddress("ttc_l1_mass",  &ev.ttc_l1_mass);
  t->SetBranchAddress("ttc_l2_pt",    &ev.ttc_l2_pt);
  t->SetBranchAddress("ttc_l2_eta",   &ev.ttc_l2_eta);
  t->SetBranchAddress("ttc_l2_phi",   &ev.ttc_l2_phi);
  t->SetBranchAddress("ttc_l2_mass",  &ev.ttc_l2_mass);
  t->SetBranchAddress("ttc_mll",      &ev.ttc_mll);
  t->SetBranchAddress("ttc_drll",     &ev.ttc_drll);
  t->SetBranchAddress("ttc_dphill",   &ev.ttc_dphill);
  t->SetBranchAddress("ttc_met",      &ev.ttc_met);
  t->SetBranchAddress("ttc_met_phi",  &ev.ttc_met_phi);
  t->SetBranchAddress("ttc_dr_l1j1",  &ev.ttc_dr_l1j1);
  t->SetBranchAddress("ttc_dr_l1j2",  &ev.ttc_dr_l1j2);
  t->SetBranchAddress("ttc_dr_l1j3",  &ev.ttc_dr_l1j3);
  t->SetBranchAddress("ttc_dr_l1j4",  &ev.ttc_dr_l1j4);
  t->SetBranchAddress("ttc_dr_l2j1",  &ev.ttc_dr_l2j1);
  t->SetBranchAddress("ttc_dr_l2j2",  &ev.ttc_dr_l2j2);
  t->SetBranchAddress("ttc_dr_l2j3",  &ev.ttc_dr_l2j3);
  t->SetBranchAddress("ttc_dr_l2j4",  &ev.ttc_dr_l2j4);
  t->SetBranchAddress("ttc_mllj1",    &ev.ttc_mllj1);
  t->SetBranchAddress("ttc_mllj2",    &ev.ttc_mllj2);
  t->SetBranchAddress("ttc_mllj3",    &ev.ttc_mllj3);
  t->SetBranchAddress("ttc_mllj4",    &ev.ttc_mllj4);

  //BASIC INFORMATION --- WZ REGION
  t->SetBranchAddress("WZ_region",    &ev.WZ_region);
  t->SetBranchAddress("WZ_zl1_id",    &ev.WZ_zl1_id);
  t->SetBranchAddress("WZ_zl2_id",    &ev.WZ_zl2_id);
  t->SetBranchAddress("WZ_wl_id",     &ev.WZ_wl_id);
  t->SetBranchAddress("WZ_zl1_pt",    &ev.WZ_zl1_pt);
  t->SetBranchAddress("WZ_zl1_eta",   &ev.WZ_zl1_eta);
  t->SetBranchAddress("WZ_zl1_phi",   &ev.WZ_zl1_phi);
  t->SetBranchAddress("WZ_zl1_mass",  &ev.WZ_zl1_mass);
  t->SetBranchAddress("WZ_zl2_pt",    &ev.WZ_zl2_pt);
  t->SetBranchAddress("WZ_zl2_eta",   &ev.WZ_zl2_eta);
  t->SetBranchAddress("WZ_zl2_phi",   &ev.WZ_zl2_phi);
  t->SetBranchAddress("WZ_zl2_mass",  &ev.WZ_zl2_mass);
  t->SetBranchAddress("WZ_l3_pt",     &ev.WZ_l3_pt);
  t->SetBranchAddress("WZ_l3_eta",    &ev.WZ_l3_eta);
  t->SetBranchAddress("WZ_l3_phi",    &ev.WZ_l3_phi);
  t->SetBranchAddress("WZ_l3_mass",   &ev.WZ_l3_mass);
  t->SetBranchAddress("WZ_Z_mass",    &ev.WZ_Z_mass);
  t->SetBranchAddress("WZ_Z_pt",      &ev.WZ_Z_pt);
  t->SetBranchAddress("WZ_Z_eta",     &ev.WZ_Z_eta);
  t->SetBranchAddress("WZ_Z_phi",     &ev.WZ_Z_phi);
  t->SetBranchAddress("WZ_met",       &ev.WZ_met);

  //BASIC INFORMATION --- DY REGION
  t->SetBranchAddress("DY_region",   &ev.DY_region);
  t->SetBranchAddress("DY_l1_id",    &ev.DY_l1_id);
  t->SetBranchAddress("DY_l2_id",    &ev.DY_l2_id);
  t->SetBranchAddress("DY_l1_pdgid", &ev.DY_l1_pdgid);
  t->SetBranchAddress("DY_l2_pdgid", &ev.DY_l2_pdgid);
  t->SetBranchAddress("DY_l1_pt",    &ev.DY_l1_pt);
  t->SetBranchAddress("DY_l1_pt_raw",&ev.DY_l1_pt_raw);
  t->SetBranchAddress("DY_l1_eta",   &ev.DY_l1_eta);
  t->SetBranchAddress("DY_l1_phi",   &ev.DY_l1_phi);
  t->SetBranchAddress("DY_l1_mass",  &ev.DY_l1_mass);
  t->SetBranchAddress("DY_l2_pt",    &ev.DY_l2_pt);
  t->SetBranchAddress("DY_l2_pt_raw",&ev.DY_l2_pt_raw);
  t->SetBranchAddress("DY_l2_eta",   &ev.DY_l2_eta);
  t->SetBranchAddress("DY_l2_phi",   &ev.DY_l2_phi);
  t->SetBranchAddress("DY_l2_mass",  &ev.DY_l2_mass);
  t->SetBranchAddress("DY_z_mass",   &ev.DY_z_mass);
  t->SetBranchAddress("DY_z_mass_raw",&ev.DY_z_mass_raw);
  t->SetBranchAddress("DY_z_pt",     &ev.DY_z_pt);
  t->SetBranchAddress("DY_z_pt_raw", &ev.DY_z_pt_raw);
  t->SetBranchAddress("DY_z_eta",    &ev.DY_z_eta);
  t->SetBranchAddress("DY_z_eta_raw",&ev.DY_z_eta_raw);
  t->SetBranchAddress("DY_z_phi",    &ev.DY_z_phi);
  t->SetBranchAddress("DY_z_phi_raw",&ev.DY_z_phi_raw);
  t->SetBranchAddress("DY_drll",     &ev.DY_drll);
}
