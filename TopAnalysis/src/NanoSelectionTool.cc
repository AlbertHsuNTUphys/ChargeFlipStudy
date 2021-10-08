#include "TopLJets2015/TopAnalysis/interface/NanoSelectionTool.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;


NanoSelectionTool::NanoSelectionTool(TTree *t, bool is_mc){

  //////////////////////
  //  INITIALIZATION  //
  //////////////////////

  t->SetBranchAddress("nMuon",                &nMuon);
  t->SetBranchAddress("Muon_tightId",         Muon_tightId);
  t->SetBranchAddress("Muon_pfRelIso04_all",  Muon_pfRelIso04_all);
  t->SetBranchAddress("Muon_eta",             Muon_eta);
  t->SetBranchAddress("Muon_tightCharge",     Muon_tightCharge);
  t->SetBranchAddress("Muon_corrected_pt",    Muon_corrected_pt);
  t->SetBranchAddress("Muon_phi",             Muon_phi);
  t->SetBranchAddress("Muon_mass",            Muon_mass);
  t->SetBranchAddress("Muon_charge",          Muon_charge);
  t->SetBranchAddress("Muon_pdgId",           Muon_pdgId);
  if(is_mc) t->SetBranchAddress("Muon_genPartIdx",      Muon_genPartIdx);

  t->SetBranchAddress("nElectron",            &nElectron);
  t->SetBranchAddress("Electron_cutBased",    Electron_cutBased);
  t->SetBranchAddress("Electron_tightCharge", Electron_tightCharge);
  t->SetBranchAddress("Electron_deltaEtaSC",  Electron_deltaEtaSC);
  t->SetBranchAddress("Electron_eta",         Electron_eta);
  t->SetBranchAddress("Electron_dxy",         Electron_dxy);
  t->SetBranchAddress("Electron_dz",          Electron_dz);
  t->SetBranchAddress("Electron_pt",          Electron_pt);
  t->SetBranchAddress("Electron_phi",         Electron_phi);
  t->SetBranchAddress("Electron_mass",        Electron_mass);
  t->SetBranchAddress("Electron_pdgId",       Electron_pdgId);
  t->SetBranchAddress("Electron_charge",      Electron_charge);
  if(is_mc) t->SetBranchAddress("Electron_genPartIdx",  Electron_genPartIdx);

  t->SetBranchAddress("nJet",                 &nJet);
  t->SetBranchAddress("Jet_pt_nom",           Jet_pt_nom);
  t->SetBranchAddress("Jet_eta",              Jet_eta);
  t->SetBranchAddress("Jet_phi",              Jet_phi);
  t->SetBranchAddress("Jet_mass_nom",         Jet_mass_nom);
  t->SetBranchAddress("Jet_jetId",            Jet_jetId);

  if(is_mc) t->SetBranchAddress("GenDressedLepton_eta",          GenDressedLepton_eta);
  if(is_mc) t->SetBranchAddress("GenDressedLepton_phi",          GenDressedLepton_phi);
  is_mc_ = is_mc;

}

std::vector<nanoParticle> NanoSelectionTool::selLeptons(){

  std::vector<nanoParticle> leptons;


  ////////////
  //  Muon  //
  ////////////
  
  for(UInt_t imu=0; imu<nMuon; imu++){

      if(!(Muon_tightId[imu])) continue;
      if( Muon_pfRelIso04_all[imu] < 0.15 && fabs(Muon_eta[imu]) < 2.4 && Muon_tightCharge[imu] == 2 &&
          Muon_corrected_pt[imu] > 15){
          Float_t unc = 0; // Uncertainty part not yet done
          TLorentzVector mp4;
          mp4.SetPtEtaPhiM(Muon_corrected_pt[imu], Muon_eta[imu], Muon_phi[imu], Muon_mass[imu]);
          
          //Gen-Matching
          int gen_idx = -1;
          if(is_mc_){
              gen_idx = Muon_genPartIdx[imu];
              if(!(gen_idx<0)){
                  if(deltaR(Muon_eta[imu],Muon_phi[imu],GenDressedLepton_eta[gen_idx],GenDressedLepton_phi[gen_idx])>0.4) gen_idx = -1;
              }
          }

          leptons.push_back(nanoParticle(mp4, Muon_charge[imu], Muon_pdgId[imu], 13, imu, gen_idx, 1.0, unc)); //puppi and qualityFlags not yet done.

      }

  }	

  ////////////////
  //  Electron  //
  ////////////////

  for(UInt_t iele=0; iele<nElectron; iele++){

      if(!(Electron_cutBased[iele]==4)) continue;
      if(Electron_tightCharge[iele]==2 && 
         ((fabs(Electron_eta[iele]+Electron_deltaEtaSC[iele]) < 1.4442 && fabs(Electron_dxy[iele]) < 0.05 && fabs(Electron_dz[iele]) < 0.1) ||
          (fabs(Electron_eta[iele]+Electron_deltaEtaSC[iele]) > 1.566 && fabs(Electron_eta[iele]+Electron_deltaEtaSC[iele]) < 2.4 && fabs(Electron_dxy[iele]) < 0.1 && fabs(Electron_dz[iele]) < 0.2)) &&
         Electron_pt[iele]>15){
          
          Float_t unc = 0.; //Uncertainty part not yet done
          TLorentzVector ep4;
          ep4.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele], Electron_phi[iele], Electron_mass[iele]);

          //Gen-Matching
          int gen_idx = -1;
          if(is_mc_){
              gen_idx = Electron_genPartIdx[iele];
              if(!(gen_idx<0)){
                  if(deltaR(Electron_eta[iele],Electron_phi[iele],GenDressedLepton_eta[gen_idx],GenDressedLepton_phi[gen_idx])>0.4) gen_idx = -1;
              }
          }

          leptons.push_back(nanoParticle(ep4, Electron_charge[iele],  Electron_pdgId[iele], 11, iele, gen_idx, 1.0, unc)); //puppi and qualityFlags not yet done.
      }

  }


  return leptons;  

}

