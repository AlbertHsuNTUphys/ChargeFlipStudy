#include "TopLJets2015/TopAnalysis/interface/NanoScaleFactorsWrapper.h"

#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>

using namespace std;

NanoScaleFactorsWrapper::NanoScaleFactorsWrapper(bool is_mc, TString era){

    if(!is_mc) return;

    if(era.Contains("era2017")) era_=2017;
    if(era.Contains("era2016")) era_=2016;

    cout << "[ScaleFactorsWrapper]" << endl;

    ////////////////////////////
    //  TRIGGER SCALE FACTOR  //
    ////////////////////////////

    TString f_triggerSF = TString(era+"/TriggerSF.root");
    TFile *fIn  = TFile::Open(f_triggerSF);
 
    cout << "Will use " << f_triggerSF << endl;

    if(fIn && !fIn->IsZombie()) {
        scaleFactorsH_["trigger_ee_l1"] = (TH2 *)fIn->Get("h2D_SF_ee_lep1pteta")->Clone();
        scaleFactorsH_["trigger_ee_l2"] = (TH2 *)fIn->Get("h2D_SF_ee_lep2pteta")->Clone();
        scaleFactorsH_["trigger_mm_l1"] = (TH2 *)fIn->Get("h2D_SF_mumu_lep1pteta")->Clone();
        scaleFactorsH_["trigger_mm_l2"] = (TH2 *)fIn->Get("h2D_SF_mumu_lep2pteta")->Clone();
        scaleFactorsH_["trigger_em_l1"] = (TH2 *)fIn->Get("h2D_SF_emu_lep1pteta")->Clone();
        scaleFactorsH_["trigger_em_l2"] = (TH2 *)fIn->Get("h2D_SF_emu_lep2pteta")->Clone();
    }    
    cout << "Complete loading Trigger scale factor" << endl;

    /////////////////////////////////////////
    //  CHARGE MISIDENTIFIED SCALE FACTOR  //
    /////////////////////////////////////////

    TString ChargeFlip_SF_chi2, ChargeFlip_SF_MLE;
    if(era_==2017){
        ChargeFlip_SF_chi2 = era+"/ChargeFlipProbability_2017_chi2_Fit.root";
        ChargeFlip_SF_MLE = era+"/ChargeFlipProbability_2017_MLE_Fit.root";
    }
    gSystem->ExpandPathName(ChargeFlip_SF_chi2);
    fIn=TFile::Open(ChargeFlip_SF_chi2);
    if(fIn && !fIn->IsZombie()) {
        cout << "electrons: charge flip SF from" << ChargeFlip_SF_chi2 << endl;
        scaleFactorsH_["ChargeFlip_chi2_data"]=(TH2F *)fIn->Get("data_CFRate");
        scaleFactorsH_["ChargeFlip_chi2_MC"]=(TH2F *)fIn->Get("MC_CFRate");
        scaleFactorsH_["ChargeFlip_chi2_data"]->SetDirectory(0);
        scaleFactorsH_["ChargeFlip_chi2_MC"]->SetDirectory(0);
    }
    fIn->Close();
    gSystem->ExpandPathName(ChargeFlip_SF_MLE);
    fIn=TFile::Open(ChargeFlip_SF_MLE);
    if(fIn && !fIn->IsZombie()) {
        cout << "electrons: charge flip SF from" << ChargeFlip_SF_MLE << endl;
        scaleFactorsH_["ChargeFlip_MLE_data"]=(TH2F *)fIn->Get("data_CFRate");
        scaleFactorsH_["ChargeFlip_MLE_MC"]=(TH2F *)fIn->Get("MC_CFRate");
        scaleFactorsH_["ChargeFlip_MLE_data"]->SetDirectory(0);
        scaleFactorsH_["ChargeFlip_MLE_MC"]->SetDirectory(0);
    }
    chargeflip_SF = 1.0;
    chargeflip_evwgt_woSF = 0.0;
    chargeflip_evwgt_SF = 0.0;
    chargeflip_count = 0;
    chargeflip_ratio_count = 0;
    fIn->Close();

}

double NanoScaleFactorsWrapper::GetEETriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta){

    if(l1_pt > 200) l1_pt = 199;
    if(l2_pt > 200) l2_pt = 199;
    
    TH2* h1_ee = scaleFactorsH_["trigger_ee_l1"];
    TH2* h2_ee = scaleFactorsH_["trigger_ee_l2"];

    Float_t sf_l1 = h1_ee->GetBinContent(h1_ee->FindBin(l1_pt,fabs(l1_eta)));
    Float_t sf_l2 = h2_ee->GetBinContent(h2_ee->FindBin(l2_pt,fabs(l2_eta)));

    return sf_l1*sf_l2;

}

double NanoScaleFactorsWrapper::GetMMTriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta){

    if(l1_pt > 200) l1_pt = 199;
    if(l2_pt > 200) l2_pt = 199;

    TH2* h1_mm = scaleFactorsH_["trigger_mm_l1"];
    TH2* h2_mm = scaleFactorsH_["trigger_mm_l2"];

    Float_t sf_l1 = h1_mm->GetBinContent(h1_mm->FindBin(l1_pt,fabs(l1_eta)));
    Float_t sf_l2 = h2_mm->GetBinContent(h2_mm->FindBin(l2_pt,fabs(l2_eta)));

    return sf_l1*sf_l2;

}

double NanoScaleFactorsWrapper::GetEMTriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta){

    if(l1_pt > 200) l1_pt = 199;
    if(l2_pt > 200) l2_pt = 199;

    TH2* h1_em = scaleFactorsH_["trigger_em_l1"];
    TH2* h2_em = scaleFactorsH_["trigger_em_l2"];

    Float_t sf_l1 = h1_em->GetBinContent(h1_em->FindBin(l1_pt,fabs(l1_eta)));
    Float_t sf_l2 = h2_em->GetBinContent(h2_em->FindBin(l2_pt,fabs(l2_eta)));

    return sf_l1*sf_l2;

}

NanoEffCorrection_t NanoScaleFactorsWrapper::GetChargeFlipScaleFactor(TString File_name, TString Fit_type, Float_t pt1,Float_t eta1, Int_t charge1, Float_t pt2, Float_t eta2, Int_t charge2){

  NanoEffCorrection_t sf(1.0,0.01); // Uncertainty part not yet be done.
  if(File_name.Contains("TTTo") || File_name.Contains("DY")){ // Only apply to ttbar and DY process
    TH2 *h_data=scaleFactorsH_["ChargeFlip_chi2_data"];
    TH2 *h_MC  =scaleFactorsH_["ChargeFlip_chi2_MC"];
    if(Fit_type.Contains("MLE")){
      TH2 *h_data_MLE = scaleFactorsH_["ChargeFlip_MLE_data"];
      TH2 *h_MC_MLE = scaleFactorsH_["ChargeFlip_MLE_MC"];
      h_data = h_data_MLE;
      h_MC = h_MC_MLE;
    }
    Float_t pt1ForSF = pt1;
    Float_t pt2ForSF = pt2;
    Float_t eta1ForSF = fabs(eta1);
    Float_t eta2ForSF = fabs(eta2);
    if (pt1 > h_data->GetXaxis()->GetXmax()) pt1ForSF = h_data->GetXaxis()->GetXmax() - 0.01;
    if (pt1 < h_data->GetXaxis()->GetXmin()) pt1ForSF = h_data->GetXaxis()->GetXmin() + 0.01;
    if (pt2 > h_data->GetXaxis()->GetXmax()) pt2ForSF = h_data->GetXaxis()->GetXmax() - 0.01;
    if (pt2 < h_data->GetXaxis()->GetXmin()) pt2ForSF = h_data->GetXaxis()->GetXmin() + 0.01;
    if (eta1 > h_data->GetYaxis()->GetXmax()) eta1ForSF = h_data->GetYaxis()->GetXmax() - 0.01;
    if (eta1 < h_data->GetYaxis()->GetXmin()) eta1ForSF = h_data->GetYaxis()->GetXmin() + 0.01;
    if (eta2 > h_data->GetYaxis()->GetXmax()) eta2ForSF = h_data->GetYaxis()->GetXmax() - 0.01;
    if (eta2 < h_data->GetYaxis()->GetXmin()) eta2ForSF = h_data->GetYaxis()->GetXmin() + 0.01;
    int   y1binForSF = h_data->GetYaxis()->FindBin(eta1ForSF);
    int   x1binForSF = h_data->GetXaxis()->FindBin(pt1ForSF);
    int   y2binForSF = h_data->GetYaxis()->FindBin(eta2ForSF);
    int   x2binForSF = h_data->GetXaxis()->FindBin(pt2ForSF);
    Float_t data_CF,MC_CF,P1_data,P2_data,P1_MC,P2_MC;
    P1_data = h_data->GetBinContent(x1binForSF,y1binForSF);
    P2_data = h_data->GetBinContent(x2binForSF,y2binForSF);
    P1_MC = h_MC->GetBinContent(x1binForSF, y1binForSF);
    P2_MC = h_MC->GetBinContent(x2binForSF, y2binForSF);
    data_CF = P1_data + P2_data - 2.*P1_data*P2_data;
    MC_CF = P1_MC + P2_MC -2.*P1_MC*P2_MC;
    if(charge1*charge2>0) sf.first = data_CF/MC_CF;
    else if(charge1*charge2<0) sf.first = (1.-data_CF)/(1.-MC_CF);
  }
  chargeflip_SF = sf.first;
  chargeflip_count++;
  return sf;
   
}


NanoScaleFactorsWrapper::~NanoScaleFactorsWrapper(){
}

