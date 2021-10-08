#ifndef _nanoselection_tool_h_
#define _nanoselection_tool_h_

#include <vector>
#include <map>
#include <iostream>

#include "TString.h"
#include "TH1.h"
#include "TTree.h"

#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

using namespace std;

class NanoSelectionTool{

  public:

  NanoSelectionTool(TTree *t, bool is_mc);
 ~NanoSelectionTool() {}
  std::vector<nanoParticle> selLeptons();


  private:

  UInt_t nMuon, nElectron, nJet;

  Bool_t Muon_tightId[500];
  Int_t  Muon_tightCharge[500], Muon_charge[500], Muon_pdgId[500], Muon_genPartIdx[500];
  Float_t  Muon_pfRelIso04_all[500], Muon_eta[500], Muon_corrected_pt[500], Muon_phi[500],
           Muon_mass[500];

  Int_t Electron_cutBased[500], Electron_tightCharge[500], Electron_charge[500], Electron_pdgId[500], Electron_genPartIdx[500];
  Float_t Electron_deltaEtaSC[500], Electron_eta[500], Electron_dxy[500], Electron_dz[500],
          Electron_pt[500], Electron_phi[500], Electron_mass[500];

  Int_t Jet_jetId[500];
  Float_t Jet_pt_nom[500], Jet_eta[500], Jet_phi[500], Jet_mass_nom[500];

  Float_t GenDressedLepton_phi[500], GenDressedLepton_eta[500];

  bool is_mc_;
};

#endif
