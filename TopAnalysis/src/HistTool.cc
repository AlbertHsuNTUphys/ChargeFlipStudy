#include "TopLJets2015/TopAnalysis/interface/HistTool.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include <iostream>

// Histogram tool for automatic creation of 2D uncertainty histograms
HistTool::HistTool(unsigned int nsyst) :
  nsyst_(nsyst)
{
  SF_need_to_norm = 1.0;
}

//
void HistTool::addHist(TString title, TH1* hist) {
  if(hist->InheritsFrom("TH2")) {
    all2dPlots_[title]=(TH2 *)hist;
    all2dPlots_[title]->SetDirectory(0);
    all2dPlots_woselSF_[title]=(TH2 *)hist;
  }
  else {
    allPlots_[title] = hist;
    allPlots_woselSF_[title] = hist;
    allPlots_[title]->SetDirectory(0);
    if (nsyst_ > 0) {
      all2dPlots_[title+"_syst"] = new TH2F(title+"_syst", hist->GetTitle(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), nsyst_+1, -0.5, nsyst_+0.5);
      all2dPlots_[title+"_syst"]->SetXTitle(hist->GetXaxis()->GetTitle());
      all2dPlots_[title+"_syst"]->SetYTitle("Variation (0=default)");
      all2dPlots_[title+"_syst"]->SetDirectory(0);
      all2dPlots_woselSF_[title+"_syst"]=(TH2 *)all2dPlots_[title+"_syst"]->Clone(title+"_syst_woSF");
    }
  }
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights,std::vector<TString> cats) {
  for(auto &c : cats)
    fill(title,value,weights,c);
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights,TString cat) {

  if (not allPlots_.count(title)) {
    //std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }

  if(allPlots_[title]->InheritsFrom("TH2")) return;

  //category specific plot, init if needed
  if(cat!=""){
    TString newTitle=cat+"_"+title;
    if(not allPlots_.count(newTitle)) {
      //std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
      allPlots_[newTitle]=(TH1 *)allPlots_[title]->Clone(newTitle);
      allPlots_[newTitle]->SetDirectory(0);
      allPlots_[newTitle]->Reset("ICE");
      allPlots_woselSF_[newTitle]=(TH1 *)allPlots_[title]->Clone(newTitle+"_woselSF");
      allPlots_woselSF_[newTitle]->Reset("ICE");
    }
    title=newTitle;
  }

  allPlots_[title]->Fill(value, weights[0]);
  allPlots_woselSF_[title]->Fill(value, weights[0]/SF_need_to_norm);

  if (nsyst_ > 0) {
    if (weights.size() > nsyst_)
      std::cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << std::endl;
    all2dPlots_[title+"_syst"]->Fill(value, 0., weights[0]);
    for (unsigned int i = 1; i < weights.size(); ++i) {
      all2dPlots_[title+"_syst"]->Fill(value, i, weights[0]*weights[i]);
    }
  }
}



void HistTool::fill2D(TString title, double valueX, double valueY, std::vector<double> weights,std::vector<TString> cats) {
  for(auto &c : cats)
    fill2D(title,valueX,valueY,weights,c);
}


void HistTool::fill2D(TString title, double valueX, double valueY, std::vector<double> weights,TString cat) {

  if (not all2dPlots_.count(title)) {
    std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }


  //category specific plot, init if needed
  if(cat!=""){
    TString newTitle=cat+"_"+title;
    if(not all2dPlots_.count(newTitle)) {
      std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
      all2dPlots_[newTitle]=(TH2 *)all2dPlots_[title]->Clone(newTitle);
      all2dPlots_[newTitle]->SetDirectory(0);
      all2dPlots_[newTitle]->Reset("ICE");
      all2dPlots_woselSF_[newTitle]=(TH2 *)all2dPlots_woselSF_[title]->Clone(newTitle);
      all2dPlots_woselSF_[newTitle]->Reset("ICE");
    }
    title=newTitle;
  }

  all2dPlots_[title]->Fill(valueX,valueY, weights[0]);
  all2dPlots_woselSF_[title]->Fill(valueX,valueY, weights[0]/SF_need_to_norm);
}

void HistTool::scale(TString title, double value,std::vector<TString> cats) {
  for(auto &c : cats)
    scale(title,value,c);
}

void HistTool::scale(TString title, double value, TString cat){
  if (not allPlots_.count(title)) {
    return;
  }
  if(allPlots_[title]->InheritsFrom("TH2")) return;
  if(cat!=""){
    TString newTitle=cat+"_"+title;
    if(not allPlots_.count(newTitle)) {
      return;
    }
    title=newTitle;
  }
  allPlots_[title]->Scale(value);
}

float HistTool::integral(TString title, TString cat){
  title = cat+"_"+title;
  if (not allPlots_.count(title)) {
    return 0.0;
  }
  return allPlots_[title]->Integral();
}

// Below are used to calculate and apply norm weight
//---------------------------------------------------------------------------

void HistTool::load_SF(double SF){
    SF_need_to_norm*=SF;
}

void HistTool::start_new_event(){
  SF_need_to_norm = 1.;
}

void HistTool::apply_norm_weight(){
    std::cout<<"--------------- 1D Histrogram -----------------"<<std::endl;
    for (std::map<TString, TH1 *>::iterator it=allPlots_.begin(); it!=allPlots_.end(); ++it){
      TString name = it->first;
      int nbins = allPlots_[name]->GetNbinsX() ;
      float total_weight = allPlots_[name]->Integral(0,nbins+1);
      float total_weight_woselSF = allPlots_woselSF_[name]->Integral(0,nbins+1);
      float scale = 1.0;
      if(total_weight!=0.0) scale = total_weight_woselSF/total_weight;
      allPlots_[it->first]->Scale(scale);
     // std::cout<<name<<" : "<<total_weight<<" | "<<total_weight_woselSF<<" | "<<allPlots_[name]->Integral(0,nbins+1)<<std::endl;
    }
   //std::cout<<"--------------- 2D Histrogram -----------------"<<std::endl;
    for (std::map<TString, TH2 *>::iterator it=all2dPlots_.begin(); it!=all2dPlots_.end(); ++it){
      TString name = it->first;
      int nxbins = all2dPlots_[name]->GetNbinsX() ;
      int nybins = all2dPlots_[name]->GetNbinsY() ;
      float total_weight = all2dPlots_[name]->Integral(0,nxbins+1,0,nybins+1);
      float total_weight_woselSF = all2dPlots_woselSF_[name]->Integral(0,nxbins+1,0,nybins+1);
      float scale = 1.0;
      if(total_weight!=0.0) scale = total_weight_woselSF/total_weight;
      all2dPlots_[it->first]->Scale(scale);
     // std::cout<<name<<" : "<<total_weight<<" | "<<total_weight_woselSF<<" | "<<all2dPlots_[name]->Integral(0,nxbins+1,0,nybins+1)<<std::endl;
    }

}
