#ifndef _hist_tool_h_
#define _hist_tool_h_

#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/CtagScaleFactorsWrapper.h"


#include <map>
#include <vector>

class HistTool {

 public:
  HistTool(unsigned int nsyst = 20);
  ~HistTool() {}

  void setNsyst(unsigned int nsyst) { nsyst_ = nsyst; }
  void addHist(TString title, TH1* hist);
  void fill(TString title, double value, std::vector<double> weights,std::vector<TString> cats);
  void fill(TString title, double value, std::vector<double> weights,TString cat="");
  void fill2D(TString title, double valueX, double valueY, std::vector<double> weights,std::vector<TString> cats);
  void fill2D(TString title, double valueX, double valueY, std::vector<double> weights,TString cat="");
  void fill(TString title, double value, double weight,TString cat="") { fill(title,value,std::vector<double>(1,weight),cat); }
  void fill(TString title, double value, double weight,std::vector<TString> cats) { fill(title,value,std::vector<double>(1,weight),cats); }

  void fill2D(TString title, double valueX, double valueY, double weight,TString cat="") { fill2D(title, valueX, valueY,std::vector<double>(1,weight),cat); }
  void fill2D(TString title, double valueX, double valueY, double weight,std::vector<TString> cats) { fill2D(title, valueX, valueY,std::vector<double>(1,weight),cats); }
  void scale(TString title, double value, TString cats);
  void scale(TString title, double value,std::vector<TString> cats);
  float integral(TString title, TString cat);

  void start_new_event();
  void apply_norm_weight();
  void load_SF(double SF);

  std::map<TString, TH1 *> &getPlots()   { return allPlots_; }
  std::map<TString, TH2 *> &get2dPlots() { return all2dPlots_; }

 private:
  unsigned int nsyst_;
  std::map<TString, TH1 *> allPlots_;
  std::map<TString, TH2 *> all2dPlots_;
  std::map<TString, TH1 *> allPlots_woselSF_;
  std::map<TString, TH2 *> all2dPlots_woselSF_;
  float SF_need_to_norm;
};

#endif
