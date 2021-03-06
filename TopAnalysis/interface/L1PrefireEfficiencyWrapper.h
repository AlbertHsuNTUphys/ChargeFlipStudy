#ifndef _L1PrefireEfficiencyWrapper_h_
#define _L1PrefireEfficiencyWrapper_h_

#include "TH2.h"
#include "TString.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"

class L1PrefireEfficiencyWrapper 
{
 public:
  L1PrefireEfficiencyWrapper(bool isData,TString era);
  EffCorrection_t getCorrection(std::vector<Jet> &jets,std::vector<Particle> &photons,bool byMax=false);
  EffCorrection_t getCorrection(std::vector<Jet> &jets,bool byMax=false);
  ~L1PrefireEfficiencyWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2 *> effMapsH_;
};

#endif
