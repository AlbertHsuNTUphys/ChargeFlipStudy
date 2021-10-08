#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/ChargeFlip_nanoaod.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/NanoEvent.h"
#include "TopLJets2015/TopAnalysis/interface/NanoSelectionTool.h"
#include "TopLJets2015/TopAnalysis/interface/NanoScaleFactorsWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "TMath.h"

using namespace std;

bool in_kinematicRegion(int i, int j, nanoParticle p, std::vector<Float_t> pt_region, std::vector<Float_t> eta_region);
bool in_kinematicRegion(int i, int j, int ii, int jj, nanoParticle p1, nanoParticle p2, std::vector<Float_t> pt_region, std::vector<Float_t> eta_region);

void RunChargeFlip_nanoaod(const TString in_fname, 
                                 TString outname,
                                 TH1F   *normH,
                                 TString era, 
                                 Bool_t  debug)
{
  //////////////////////
  //  INITIALIZATION  //
  //////////////////////

  TString dirin  = gSystem->DirName(in_fname);
  TString baseMC = gSystem->BaseName(dirin); cout<<"Base MC name: "<<baseMC<<endl;

  //PREPARE OUTPUT FILE

  TString baseName = gSystem->BaseName(outname);
  TString dirName  = gSystem->DirName(outname);
  TFile *fOut      = TFile::Open(dirName+"/"+baseName, "RECREATE");

  fOut -> cd();

  //CREATE OUPUT TREE

  const char* treename = "TreeInput";
  TTree t_output(treename, treename);

  //INPUTFILE
 
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()){
      cout << "Corrupted or missing file " << in_fname << endl;
      return;
  }
  bool is_mc = outname.Contains("MC");

  //INPUTTREE  
  
  TTree *t = (TTree*)f->Get("Events");
  
  NanoEvent_t ev;
  attachToNanoEventTree(t,ev,true,is_mc);

  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(1000,nentries);

  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;


  //////////////////////
  //  ANALYSIS SETUP  //
  //////////////////////

  NanoSelectionTool selector(t, is_mc);

  NanoScaleFactorsWrapper SFWrapper(is_mc, era);
  
  bool  use_nanoaod_region = 1;
  Float_t normWeight = 1.;

  if(is_mc){

     TH1F *NormH = (TH1F*) f->Get("nEventsGenWeighted");
     normWeight    = 1./(NormH->GetBinContent(1));

  }

  /////////////////////////
  //  CHARGE FLIP REGION //
  /////////////////////////
  
  std::vector<Float_t> pt_region  = {20., 50.,  100., 100000000000.};
  std::vector<Float_t> eta_region = {0.,  0.8, 1.479, 2.4};
  int pt_bins  = pt_region.size() - 1;
  int eta_bins = eta_region.size() -1 ;

  Float_t OS_Zmass = 90.44464129969609;
  Float_t OS_width = 13.24949999998246;
  Float_t SS_Zmass = 88.92155676987167;
  Float_t SS_width = 15.675899999976805;


  Float_t t_Nss[pt_bins][eta_bins][pt_bins][eta_bins], t_Nos[pt_bins][eta_bins][pt_bins][eta_bins];
  Float_t t_g_Nss[pt_bins][eta_bins], t_g_Nos[pt_bins][eta_bins];
  Float_t t_gm_Nss[pt_bins][eta_bins], t_gm_Nos[pt_bins][eta_bins];

  for(int i = 0; i<pt_bins;i++){
    for(int j = 0; j<eta_bins; j++){
      for(int ii = 0; ii<pt_bins;ii++){
        for(int jj = 0; jj<eta_bins;jj++){
         t_Nss[i][j][ii][jj] = 0.;
         t_Nos[i][j][ii][jj] = 0.;
  }}}}
  for(int i = 0; i < pt_bins; i++){
    for(int j = 0; j < eta_bins; j++){
      t_g_Nss[i][j] = 0.;
      t_gm_Nss[i][j] = 0.;
      t_g_Nos[i][j] = 0.;
      t_gm_Nos[i][j] = 0.;
  }}




  //////////////////////
  //  BOOK HISTOGRAM  //
  //////////////////////

  fOut -> cd();
  HistTool ht;
  ht.setNsyst(0);

  // CONTROL PLOT
  ht.addHist("control_genweight", new TH1F("control_genweight",";weight      ;Events(Raw)", 20,-2,2));
  ht.addHist("control_pdgid_check", new TH1F("control_pdgid_check", ";pdgid  ;Events",    60,-30,30));
  ht.addHist("control_genmatching", new TH1F("control_genmatching", ";Matched ;Events",   4,-2,2));
  ht.addHist("control_l1_pdgid",    new TH1F("control_l1_pdgid",    ";pdgid   ;Events",   15,0,15));
  ht.addHist("control_l2_pdgid",    new TH1F("control_l2_pdgid",    ";pdgid   ;Events",   15,0,15));
  ht.addHist("control_ll_pdgid",    new TH1F("control_ll_pdgid",    ";pdgid   ;Events",   30,0,30));

  ht.addHist("check_l1_pt",       new TH1F("check_l1_pt",      ";Pt(l1)      ;Events", 60,0,60));
  ht.addHist("check_l2_pt",       new TH1F("check_l2_pt",      ";Pt(l2)      ;Events", 60,0,60));
  ht.addHist("check_nleptons",    new TH1F("check_nleptons",   ";n(lep)      ;Events", 6,-3, 3));

  ht.addHist("DY_region",         new TH1F("DY_region",        ";Region      ;Events", 5,0,5));
  ht.addHist("DY_z_mass",         new TH1F("DY_z_mass",        ";M(Z) [GeV]  ;Events", 60, 60, 120));
  ht.addHist("DY_z_pt",           new TH1F("DY_z_pt",          ";Pt(Z) [GeV] ;Events", 30,  0, 600));
  ht.addHist("DY_l1_pt",          new TH1F("DY_l1_pt",         ";Pt(l1) [GeV];Events", 30,  0, 600));
  ht.addHist("DY_l2_pt",          new TH1F("DY_l2_pt",         ";Pt(l2) [GeV];Events", 30,  0, 600));
  ht.addHist("DY_l1_eta",         new TH1F("DY_l1_eta",        ";#eta(l1)    ;Events", 10,-2.5,2.5));
  ht.addHist("DY_l2_eta",         new TH1F("DY_l2_eta",        ";#eta(l2)    ;Events", 10,-2.5,2.5));
  ht.addHist("DY_l1_pt_test",     new TH1F("DY_l1_pt_test",    ";Pt(l1) [GeV];Events", 30,  0, 600));
  ht.addHist("DY_ll_pdgId",       new TH1F("DY_ll_pdgId",      ";pdgId;       Events", 30,  0,  30));

  // CHARGE FLIP STUDY
  ht.addHist("os_Zmass", new TH1F("os_Zmass", ";M(Z) [GeV]; Events", 80, 50, 130));
  ht.addHist("os_l1_pt", new TH1F("os_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("os_l2_pt", new TH1F("os_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("os_l1_eta",new TH1F("os_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("os_l2_eta",new TH1F("os_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("os_njets", new TH1F("os_njets", ";njets;  Events", 10, 0, 10));
  ht.addHist("ss_Zmass", new TH1F("ss_Zmass", ";M(Z) [GeV]; Events", 80, 50, 130));
  ht.addHist("ss_l1_pt", new TH1F("ss_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("ss_l2_pt", new TH1F("ss_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("ss_l1_eta",new TH1F("ss_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("ss_l2_eta",new TH1F("ss_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("ss_njets", new TH1F("ss_njets", ";njets;  Events", 10, 0, 10));

  ht.addHist("os_SF_Zmass", new TH1F("os_SF_Zmass", ";M(Z) [GeV]; Events", 80, 50, 130));
  ht.addHist("os_SF_l1_pt", new TH1F("os_SF_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("os_SF_l2_pt", new TH1F("os_SF_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("os_SF_l1_eta",new TH1F("os_SF_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("os_SF_l2_eta",new TH1F("os_SF_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("os_SF_njets", new TH1F("os_SF_njets", ";njets;  Events", 10, 0, 10));
  ht.addHist("ss_SF_Zmass", new TH1F("ss_SF_Zmass", ";M(Z) [GeV]; Events", 80, 50, 130));
  ht.addHist("ss_SF_l1_pt", new TH1F("ss_SF_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("ss_SF_l2_pt", new TH1F("ss_SF_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("ss_SF_l1_eta",new TH1F("ss_SF_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("ss_SF_l2_eta",new TH1F("ss_SF_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("ss_SF_njets", new TH1F("ss_SF_njets", ";njets;  Events", 10, 0, 10));

  ht.addHist("SB_os_l1_pt", new TH1F("SB_os_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_os_l2_pt", new TH1F("SB_os_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_os_l1_eta",new TH1F("SB_os_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_os_l2_eta",new TH1F("SB_os_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_os_njets", new TH1F("SB_os_njets", ";njets;  Events", 10, 0, 10));
  ht.addHist("SB_ss_l1_pt", new TH1F("SB_ss_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_ss_l2_pt", new TH1F("SB_ss_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_ss_l1_eta",new TH1F("SB_ss_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_ss_l2_eta",new TH1F("SB_ss_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_ss_njets", new TH1F("SB_ss_njets", ";njets;  Events", 10, 0, 10));

  ht.addHist("SB_os_SF_l1_pt", new TH1F("SB_os_SF_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_os_SF_l2_pt", new TH1F("SB_os_SF_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_os_SF_l1_eta",new TH1F("SB_os_SF_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_os_SF_l2_eta",new TH1F("SB_os_SF_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_os_SF_njets", new TH1F("SB_os_SF_njets", ";njets;  Events", 10, 0, 10));
  ht.addHist("SB_ss_SF_l1_pt", new TH1F("SB_ss_SF_l1_pt", ";Pt(l1) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_ss_SF_l2_pt", new TH1F("SB_ss_SF_l2_pt", ";Pt(l2) [GeV]; Events",30, 0, 600));
  ht.addHist("SB_ss_SF_l1_eta",new TH1F("SB_ss_SF_l1_eta",";#eta(l1) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_ss_SF_l2_eta",new TH1F("SB_ss_SF_l2_eta",";#eta(l2) ; Events", 10,-2.5, 2.5));
  ht.addHist("SB_ss_SF_njets", new TH1F("SB_ss_SF_njets", ";njets;  Events", 10, 0, 10));


  /////////////////////////
  //  BOOK OUPUT BRANCH  //
  /////////////////////////

  int event_index;
  Int_t  l1_pid, l2_pid;
  Int_t l1_id,  l2_id;
  Int_t l3_pid;
  Int_t l3_id = -1;

  t_output.Branch("DY_z_mass",   &ev.DY_z_mass, "DY_Z_mass/F");
  t_output.Branch("event_index", &event_index,  "event_index/I");
  t_output.Branch("l1_pid",      &l1_pid,       "l1_pid/I");
  t_output.Branch("l2_pid",      &l2_pid,       "l2_pid/I");
  t_output.Branch("l3_pid",      &l3_pid,       "l3_pid/I");
  t_output.Branch("l1_id",       &l1_id,        "l1_id/I");
  t_output.Branch("l2_id",       &l2_id,        "l2_id/I");
  t_output.Branch("l3_id",       &l3_id,        "l3_id/I");


  //////////////////
  //  EVENT LOOP  //
  //////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  //  PRESELCTIONS IN NTUPLE PRODUCTION: (in 2021.09.08 version)                       //
  //  1. veto PV_npvsGood < 1                                                          //
  //  2. at least 2 tight leptons                                                      //
  //  3. satisfy at least one of following requirements:                               //
  //     i. exact 2 tight leptons with pt > 20 GeV without any loose lepton (ttc_nl).  //
  //    ii. only  3 tight leptons with pt > 20 GeV without any loose lepton,           //
  //        no b-jet, mll>4, |Z-91.1876|<15 (WZ_region)                                //
  //   iii. exact 2 opposite sign eles/muons with pt > 20 GeV. veto 3rd lepton,        //
  //        with |mll-91.1876|<15                                                      //
  ///////////////////////////////////////////////////////////////////////////////////////

  for (Int_t iev=0; iev<nentries; iev++){
 
      t->GetEntry(iev);
      event_index = iev;
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(Float_t)iev/(Float_t)nentries); fflush(stdout); }
      ht.start_new_event();

      ///////////////
      //  TRIGGER  //
      ///////////////

      bool passtrigger_ee  = 0;
      bool passtrigger_mm  = 0;
      bool passtrigger_emu = 0;

      if(era.Contains("2017")){

          int a = ev.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
          int b = ev.HLT_IsoMu27;
          int c = ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
          int d = ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
          int e = ev.HLT_passEle32WPTight;
          int f = ev.HLT_Ele35_WPTight_Gsf;
          int g = ev.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
          int h = ev.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
          int i = ev.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

          if(outname.Contains("DoubleMuon",     TString::kIgnoreCase)){
              if(a) passtrigger_mm = 1;
          }

          if(outname.Contains("SingleMuon",     TString::kIgnoreCase)){
              if(!(a) && b)       passtrigger_mm  = 1;
              if(!(g||h||i) && b) passtrigger_emu = 1;
          }

          if(outname.Contains("DoubleEG",       TString::kIgnoreCase)){
              if(c||d) passtrigger_ee = 1;
          }

          if(outname.Contains("SingleElectron", TString::kIgnoreCase)){
              if(!(c||d) && (e||f))    passtrigger_ee  = 1;
              if(!(g||h||i) && (e||f)) passtrigger_emu = 1;
          }

          if(outname.Contains("MuonEG",         TString::kIgnoreCase)){
              if(g||h||i) passtrigger_emu = 1;
          }

          if(outname.Contains("MC",             TString::kIgnoreCase)){
              if(a||b||c||d||e||f||g||h||i){
                  passtrigger_ee  = 1;
                  passtrigger_mm  = 1;
                  passtrigger_emu = 1;
              }
          }

      }

      //////////////
      //  FILTER  //
      //////////////

      if (!(ev.Flag_HBHENoiseFilter &&
            ev.Flag_HBHENoiseIsoFilter &&
            ev.Flag_EcalDeadCellTriggerPrimitiveFilter &&
            ev.Flag_goodVertices &&
            ev.Flag_globalSuperTightHalo2016Filter &&
            ev.Flag_BadPFMuonFilter)) continue;

      /////////////////////
      //  SELECT OBJECT  //
      /////////////////////

      std::vector<nanoParticle> leptons;

      if(!use_nanoaod_region)   leptons = selector.selLeptons();
      else {

          if(!(ev.ttc_region > 0 || ev.DY_region > 0)) continue;
          TLorentzVector p4;
          if(ev.ttc_region > 0){
              p4.SetPtEtaPhiM(ev.ttc_l1_pt, ev.ttc_l1_eta, ev.ttc_l1_phi, ev.ttc_l1_mass);
              leptons.push_back(nanoParticle(p4, ev.ttc_l1_pdgid/fabs(ev.ttc_l1_pdgid), ev.ttc_l1_pdgid,
                                             fabs(ev.ttc_l1_pdgid), ev.ttc_l1_id, -1, 1.0, 0.));
              p4.SetPtEtaPhiM(ev.ttc_l2_pt, ev.ttc_l2_eta, ev.ttc_l2_phi, ev.ttc_l2_mass);
              leptons.push_back(nanoParticle(p4, ev.ttc_l2_pdgid/fabs(ev.ttc_l2_pdgid), ev.ttc_l2_pdgid, 
                                             fabs(ev.ttc_l2_pdgid), ev.ttc_l2_id, -1, 1.0, 0.));
          }

          else{
              p4.SetPtEtaPhiM(ev.DY_l1_pt, ev.DY_l1_eta, ev.DY_l1_phi, ev.DY_l1_mass);
              leptons.push_back(nanoParticle(p4, ev.DY_l1_pdgid/fabs(ev.DY_l1_pdgid), ev.DY_l1_pdgid, 
                                             fabs(ev.DY_l1_pdgid), ev.DY_l1_id, -1, 1.0, 0.));
              p4.SetPtEtaPhiM(ev.DY_l2_pt, ev.DY_l2_eta, ev.DY_l2_phi, ev.DY_l2_mass);
              leptons.push_back(nanoParticle(p4, ev.DY_l2_pdgid/fabs(ev.DY_l2_pdgid), ev.DY_l2_pdgid, 
                                             fabs(ev.DY_l2_pdgid), ev.DY_l2_id, -1, 1.0, 0.));
          }
      }

/*      sort(leptons.begin(),leptons.end(),
        [](const nanoParticle& a, const nanoParticle& b){
        return a.Pt() > b.Pt();
        }
      );*/


      if(leptons.size()<2) continue;

      ///////////////
      //  CHANNEL  //
      ///////////////

      int dimuon_event     = 0;
      int dielectron_event = 0;
      int emu_event        = 0;
      int mue_event        = 0;

      if (outname.Contains("DoubleMuon",TString::kIgnoreCase)){
         if (!((fabs(leptons[0].id()) == 13 && fabs(leptons[1].id()) == 13) && passtrigger_mm)) continue;
            dimuon_event=1;
      }

      if (outname.Contains("DoubleEG",TString::kIgnoreCase)){
          if (!((fabs(leptons[0].id()) == 11 && fabs(leptons[1].id()) == 11) && passtrigger_ee)) continue;
          dielectron_event = 1;
      }

      if (outname.Contains("MuonEG",TString::kIgnoreCase)){
         if (!((fabs(leptons[0].id())+fabs(leptons[1].id()))==24 && passtrigger_emu)) continue;
         if (fabs(leptons[0].id()) == 11 && fabs(leptons[1].id()) == 13) emu_event = 1;
         if (fabs(leptons[0].id()) == 13 && fabs(leptons[1].id()) == 11) mue_event = 1;
      }

      if (outname.Contains("SingleElectron",TString::kIgnoreCase)){
         if (!(((fabs(leptons[0].id()) == 11 && fabs(leptons[1].id()) == 11) && passtrigger_ee) ||
               ((fabs(leptons[0].id())+fabs(leptons[1].id()))==24 && passtrigger_emu))) continue;
         if ((fabs(leptons[0].id())+fabs(leptons[1].id()))==24) emu_event = 1;
         else dielectron_event = 1;
      }

      if (outname.Contains("SingleMuon",TString::kIgnoreCase)){
         if (!(((fabs(leptons[0].id()) == 13 && fabs(leptons[1].id()) == 13) && passtrigger_mm) ||
               ((fabs(leptons[0].id())+fabs(leptons[1].id()))==24 && passtrigger_emu))) continue;
         if ((fabs(leptons[0].id())+fabs(leptons[1].id()))==24) emu_event = 1;
         else dimuon_event = 1;
      }

      if (outname.Contains("MC",TString::kIgnoreCase)){
         if ((fabs(leptons[0].id()) == 11 && fabs(leptons[1].id()) == 11) && passtrigger_ee) dielectron_event = 1;
         if ((fabs(leptons[0].id()) == 13 && fabs(leptons[1].id()) == 13) && passtrigger_mm) dimuon_event = 1;
         if (((fabs(leptons[0].id()) + fabs(leptons[1].id())) == 24) && passtrigger_emu) emu_event = 1;
      }

      std::vector<TString> tags2={"inc"};
      if (dimuon_event)           tags2.push_back("mm");
      if (dielectron_event)       tags2.push_back("ee");
      if (emu_event || mue_event) tags2.push_back("emu");

      //////////////
      //  WEIGHT  //  
      //////////////


      /////////////////////////////////////////////////////////////////
      //  CORRECTION ALREADY APPLIED IN NTUPLE (version 2021.09.08)  //
      //  1. EGAMMA CORRECTION (Electron)                            //
      //  2. Rochester CORRECTION  (Muon)                            //
      /////////////////////////////////////////////////////////////////
        
      Float_t evWgt = 1.0;
      Float_t ChargeFlip_SF = 1.0;

      if (outname.Contains("MC",TString::kIgnoreCase)){

          Float_t SF_FirstLepton  = 1.;
          Float_t SF_SecondLepton = 1.;
          Float_t PrefireWeight   = 1.;
          Float_t puWeight        = 1.;
          double trigger_SF      = 1.;

          // ID SCALE FACTOR

          int idx_l1 = leptons[0].originalReference();
          int idx_l2 = leptons[1].originalReference();

          if(fabs(leptons[0].id())==11){
              SF_FirstLepton*=ev.Electron_CutBased_TightID_SF[idx_l1];
              SF_FirstLepton*=ev.Electron_RECO_SF[idx_l1];
          }
          else if(fabs(leptons[0].id())==13){
              SF_FirstLepton*=ev.Muon_CutBased_TightID_SF[idx_l1];
              SF_FirstLepton*=ev.Muon_TightISO_SF[idx_l1];
          }

          if(fabs(leptons[1].id())==11){
              SF_SecondLepton*=ev.Electron_CutBased_TightID_SF[idx_l2];
              SF_SecondLepton*=ev.Electron_RECO_SF[idx_l2];
          }
          else if(fabs(leptons[1].id())==13){
              SF_SecondLepton*=ev.Muon_CutBased_TightID_SF[idx_l2];
              SF_SecondLepton*=ev.Muon_TightISO_SF[idx_l2];
          }
        
          // L1-PREFIRING SCALE FACTOR
          
          PrefireWeight = ev.PrefireWeight;

          // PILE UP SCALE FACTOR

          puWeight      = ev.puWeight;

          // TRIGGER SCALE FACTOR

          if(dielectron_event)     trigger_SF *= SFWrapper.GetEETriggerScaleFactor(leptons[0].pt(), leptons[1].pt(), leptons[0].eta(), leptons[1].eta());
          if(dimuon_event)         trigger_SF *= SFWrapper.GetMMTriggerScaleFactor(leptons[0].pt(), leptons[1].pt(), leptons[0].eta(), leptons[1].eta());
          if(emu_event||mue_event) trigger_SF *= SFWrapper.GetEMTriggerScaleFactor(leptons[0].pt(), leptons[1].pt(), leptons[0].eta(), leptons[1].eta());

          // CHARGE FLIP SCALE FACTOR

          if(dielectron_event)
             ChargeFlip_SF *= SFWrapper.GetChargeFlipScaleFactor(in_fname,"Chi2",leptons[0].pt(),leptons[0].eta(),leptons[0].charge(),leptons[1].pt(),leptons[1].eta(),leptons[1].charge()).first;

          // WEIGHT

          evWgt*=ev.genWeight*SF_FirstLepton*SF_SecondLepton*PrefireWeight*puWeight*normWeight*trigger_SF/fabs(ev.genWeight);
      }

      //////////////////////////////////
      //  STAGE0 - CONTROL DY REGION  //    
      //////////////////////////////////

      ht.fill("control_pdgid_check", fabs(leptons[0].id())-leptons[0].nanoid(), evWgt, tags2);
      ht.fill("control_pdgid_check", fabs(leptons[1].id())-leptons[1].nanoid(), evWgt, tags2);
      ht.fill("control_l1_pdgid", fabs(leptons[0].id()), evWgt, tags2);
      ht.fill("control_l2_pdgid", fabs(leptons[1].id()), evWgt, tags2);
      ht.fill("control_ll_pdgid", fabs(leptons[0].id())+fabs(leptons[1].id()), evWgt, tags2);
      ht.fill("DY_region", ev.DY_region, evWgt, tags2);

      if (is_mc){
          ht.fill("control_genweight", ev.genWeight,1., tags2);
          if(leptons[0].genIdx()<0) ht.fill("control_genmatching", -1, evWgt, tags2);
          else                      ht.fill("control_genmatching",  1, evWgt, tags2);
          if(leptons[1].genIdx()<0) ht.fill("control_genmatching", -1, evWgt, tags2);
          else                      ht.fill("control_genmatching",  1, evWgt, tags2);
      }

      if(!(ev.ttc_region > 0)){
          ht.fill("check_l1_pt", fabs(leptons[0].pt() - ev.DY_l1_pt), evWgt, tags2);
          ht.fill("check_l2_pt", fabs(leptons[1].pt() - ev.DY_l2_pt), evWgt, tags2);
          ht.fill("check_nleptons", (leptons.size() - ev.n_tight_ele - ev.n_tight_muon), evWgt, tags2);
      }

      bool DY_Filter = ev.DY_z_mass>60. && ev.DY_z_mass<120. && (ev.DY_l1_pt>30. || ev.DY_l2_pt>30.) && ev.DY_drll>0.3;
      if ((ev.DY_region > 0 && DY_Filter)) {
          ht.fill("DY_z_mass", ev.DY_z_mass, evWgt, tags2);
          ht.fill("DY_z_pt",   ev.DY_z_pt,   evWgt, tags2);
          ht.fill("DY_l1_pt",  ev.DY_l1_pt,  evWgt, tags2);
          ht.fill("DY_l2_pt",  ev.DY_l2_pt,  evWgt, tags2);
          ht.fill("DY_l1_eta", ev.DY_l1_eta, evWgt, tags2);
          ht.fill("DY_l2_eta", ev.DY_l2_eta, evWgt, tags2);
          ht.fill("DY_l1_pt_test", leptons[0].pt(), evWgt, tags2);
      }
      if (ev.DY_region == 1) ht.fill("DY_ll_pdgId", fabs(ev.DY_l1_pdgid) + fabs(ev.DY_l2_pdgid), evWgt, tags2);

      ////////////////////////////////
      //  STAGE1 - BASIC SELECTION  //
      ////////////////////////////////

      if (!((ev.n_tight_ele + ev.n_tight_muon)==2)) continue;
      if ((ev.n_loose_ele + ev.n_loose_muon) >0)    continue;
      if (ev.n_tight_jet > 2 )                      continue; 
      if (ev.MET_pt      > 50)                      continue;
      
      // Z MASS HISTOGRAM IN OS/SS CHANNEL
      Float_t zmass = (leptons[0]+leptons[1]).M();
      if (leptons[0].charge()*leptons[1].charge() < 0) ht.fill("os_Zmass", zmass, evWgt, tags2);
      else                                             ht.fill("ss_Zmass", zmass, evWgt, tags2);


      //////////////////////////////////
      //  STAGE2 - CHARGE FLIP STUDY  //
      //////////////////////////////////   

      if (dielectron_event){
          for(int i = 0; i<pt_bins; i++){
           for(int j = 0; j<eta_bins; j++){
            for(int ii = 0; ii<pt_bins; ii++){
             for(int jj = 0; jj<eta_bins; jj++){
                 if(in_kinematicRegion(i,j,ii,jj,leptons[0],leptons[1],pt_region,eta_region)){
                     if(leptons[0].charge()*leptons[1].charge() > 0){
                         if(fabs(zmass-SS_Zmass)<SS_width) t_Nss[i][j][ii][jj]+=evWgt;             
                     }
                     else if(leptons[0].charge()*leptons[1].charge() < 0){
                           if(fabs(zmass-OS_Zmass)<OS_width) t_Nos[i][j][ii][jj]+=evWgt;
                     }
                 }
          }}}}
      }

      if (is_mc && 
          ((leptons[0].charge()*leptons[1].charge()<0 && fabs(zmass-OS_Zmass)<OS_width)||
            (leptons[0].charge()*leptons[1].charge()>0 && fabs(zmass-SS_Zmass)<SS_width))){
          for(int i = 0; i<2; i++){

            int gen_idx = leptons[i].genIdx();
            if(gen_idx<0) continue;

            for(int ii=0; ii<pt_bins; ii++){
             for(int jj=0; jj<eta_bins; jj++){
                if(!(in_kinematicRegion(ii,jj,leptons[i],pt_region,eta_region))) continue;
                if(fabs(leptons[i].id())==11){
                    if(leptons[i].id()*ev.GenDressedLepton_pdgId[gen_idx] > 0) t_g_Nss[ii][jj]+=evWgt;
                    else if (leptons[i].id()*ev.GenDressedLepton_pdgId[gen_idx] < 0) t_g_Nos[ii][jj]+=evWgt;
                }
                else if(fabs(leptons[i].id())==13){
                    if(leptons[i].id()*ev.GenDressedLepton_pdgId[gen_idx] > 0) t_gm_Nss[ii][jj]+=evWgt;
                    else if (leptons[i].id()*ev.GenDressedLepton_pdgId[gen_idx] < 0) t_gm_Nos[ii][jj]+=evWgt;
                }
            }}
          }
      }

      // TEST CHARGE FLIP SCALE FACTOR
    
      if(leptons[0].charge()*leptons[1].charge() < 0){
          ht.fill("os_l1_pt", leptons[0].pt(), evWgt, tags2);
          ht.fill("os_l2_pt", leptons[1].pt(), evWgt, tags2);
          ht.fill("os_l1_eta",leptons[0].eta(),evWgt, tags2);
          ht.fill("os_l2_eta",leptons[1].eta(),evWgt, tags2);
          ht.fill("os_njets", ev.n_tight_jet,  evWgt, tags2);
          if(!(fabs(zmass-OS_Zmass)<OS_width)){
              ht.fill("SB_os_l1_pt", leptons[0].pt(), evWgt, tags2);
              ht.fill("SB_os_l2_pt", leptons[1].pt(), evWgt, tags2);
              ht.fill("SB_os_l1_eta",leptons[0].eta(),evWgt, tags2);
              ht.fill("SB_os_l2_eta",leptons[1].eta(),evWgt, tags2);
              ht.fill("SB_os_njets", ev.n_tight_jet,  evWgt, tags2);
          }
      }

      if(leptons[0].charge()*leptons[1].charge() > 0){
          ht.fill("ss_l1_pt", leptons[0].pt(), evWgt, tags2);
          ht.fill("ss_l2_pt", leptons[1].pt(), evWgt, tags2);
          ht.fill("ss_l1_eta",leptons[0].eta(),evWgt, tags2);
          ht.fill("ss_l2_eta",leptons[1].eta(),evWgt, tags2);
          ht.fill("ss_njets", ev.n_tight_jet,  evWgt, tags2);
          if(!(fabs(zmass-SS_Zmass)<SS_width)){
              ht.fill("SB_ss_l1_pt", leptons[0].pt(), evWgt, tags2);
              ht.fill("SB_ss_l2_pt", leptons[1].pt(), evWgt, tags2);
              ht.fill("SB_ss_l1_eta",leptons[0].eta(),evWgt, tags2);
              ht.fill("SB_ss_l2_eta",leptons[1].eta(),evWgt, tags2);
              ht.fill("SB_ss_njets", ev.n_tight_jet,  evWgt, tags2);
          }
      }

      evWgt *= ChargeFlip_SF;
      ht.load_SF(ChargeFlip_SF);

      if(leptons[0].charge()*leptons[1].charge() < 0){
          ht.fill("os_SF_Zmass", zmass, evWgt, tags2);
          ht.fill("os_SF_l1_pt", leptons[0].pt(), evWgt, tags2);
          ht.fill("os_SF_l2_pt", leptons[1].pt(), evWgt, tags2);
          ht.fill("os_SF_l1_eta",leptons[0].eta(),evWgt, tags2);
          ht.fill("os_SF_l2_eta",leptons[1].eta(),evWgt, tags2);
          ht.fill("os_SF_njets", ev.n_tight_jet,  evWgt, tags2);
          if(!(fabs(zmass-OS_Zmass)<OS_width)){
              ht.fill("SB_os_SF_l1_pt", leptons[0].pt(), evWgt, tags2);
              ht.fill("SB_os_SF_l2_pt", leptons[1].pt(), evWgt, tags2);
              ht.fill("SB_os_SF_l1_eta",leptons[0].eta(),evWgt, tags2);
              ht.fill("SB_os_SF_l2_eta",leptons[1].eta(),evWgt, tags2);
              ht.fill("SB_os_SF_njets", ev.n_tight_jet,  evWgt, tags2);
          }
      }

      if(leptons[0].charge()*leptons[1].charge() > 0){
          ht.fill("ss_SF_Zmass", zmass, evWgt, tags2);
          ht.fill("ss_SF_l1_pt", leptons[0].pt(), evWgt, tags2);
          ht.fill("ss_SF_l2_pt", leptons[1].pt(), evWgt, tags2);
          ht.fill("ss_SF_l1_eta",leptons[0].eta(),evWgt, tags2);
          ht.fill("ss_SF_l2_eta",leptons[1].eta(),evWgt, tags2);
          ht.fill("ss_SF_njets", ev.n_tight_jet,  evWgt, tags2);
          if(!(fabs(zmass-SS_Zmass)<SS_width)){
              ht.fill("SB_ss_SF_l1_pt", leptons[0].pt(), evWgt, tags2);
              ht.fill("SB_ss_SF_l2_pt", leptons[1].pt(), evWgt, tags2);
              ht.fill("SB_ss_SF_l1_eta",leptons[0].eta(),evWgt, tags2);
              ht.fill("SB_ss_SF_l2_eta",leptons[1].eta(),evWgt, tags2);
              ht.fill("SB_ss_SF_njets", ev.n_tight_jet,  evWgt, tags2);
          }
      }

      evWgt /= ChargeFlip_SF;
      ht.load_SF(1./ChargeFlip_SF);

      if((leptons.size() - ev.n_tight_muon - ev.n_tight_ele)==0) continue;
      l1_pid = leptons[0].id();
      l2_pid = leptons[1].id();
      l1_id  = leptons[0].originalReference();
      l2_id  = leptons[1].originalReference();
      l3_id  = -1;
      if(leptons.size()>2){
          l3_id = leptons[2].originalReference();
          l3_pid = leptons[2].id();
      }

      /////////////////
      //  END STAGE  //
      /////////////////
      
      t_output.Fill();

  }

      ///////////////
      //  STORAGE  //
      ///////////////

      cout << "Closing input file" << endl;
      f->Close();

      cout << "Writing..." << endl;

      fOut->cd();
     
      cout << "cd fOut" << endl;

      t_output.Fill();

      TString s4 = TString("t_Nss["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");
      TString s5 = TString("t_Nos["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");
      TString s6 = TString("t_g_Nss["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");
      TString s7 = TString("t_g_Noc["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");
      TString s8 = TString("t_gm_Nss["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");
      TString s9 = TString("t_gm_Noc["+std::to_string(pt_bins)+"]["+std::to_string(eta_bins)+"]/F");

      TBranch* b_ss = t_output.Branch("t_Nss", t_Nss, s4);
      TBranch* b_os = t_output.Branch("t_Nos", t_Nos, s5);
      TBranch* b_g_ss = t_output.Branch("t_g_Nss", t_g_Nss, s6);
      TBranch* b_g_os = t_output.Branch("t_g_Nos", t_g_Nos, s7);
      TBranch* b_gm_ss = t_output.Branch("t_gm_Nss", t_gm_Nss, s8);
      TBranch* b_gm_os = t_output.Branch("t_gm_Nos", t_gm_Nos, s9);

      b_ss->Fill();
      b_os->Fill();
      b_g_ss->Fill();
      b_g_os->Fill();
      b_gm_ss->Fill();
      b_gm_os->Fill();

      t_output.Write();
   
      cout << "Add hsitogram" << endl;
      ht.apply_norm_weight();

      for (auto& it : ht.getPlots()){
          cout << "Entry" << endl;
          if(it.second->GetEntries()==0) continue;
          cout << "saving histogram " << endl;
          it.second->SetDirectory(fOut); it.second->Write();
      }
      for (auto& it : ht.get2dPlots()){
          if(it.second->GetEntries()==0) continue;
          it.second->SetDirectory(fOut); it.second->Write();
      }
   
     cout << "Saving..." << endl;
      
     fOut->Close();
}

bool in_kinematicRegion(int i, int j, nanoParticle p, std::vector<Float_t> pt_region, std::vector<Float_t> eta_region){

    bool pass = true; 
    if(!(pt_region[i] <= p.pt() && p.pt() < pt_region[i+1])) pass = false;
    else if(!(eta_region[j] <= fabs(p.eta()) && fabs(p.eta())<eta_region[j+1])) pass = false;
    return pass;

}
bool in_kinematicRegion(int i, int j, int ii, int jj, nanoParticle p1, nanoParticle p2, std::vector<Float_t> pt_region, std::vector<Float_t> eta_region){

    return in_kinematicRegion(i,j,p1,pt_region,eta_region) && in_kinematicRegion(ii,jj,p2,pt_region,eta_region);

}

