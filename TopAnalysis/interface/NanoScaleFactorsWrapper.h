#ifndef _NanoScaleFactorsWrapper_h_
#define _NanoScaleFactorsWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

typedef std::pair<double,double> NanoEffCorrection_t;

class NanoScaleFactorsWrapper
{
 public:
  NanoScaleFactorsWrapper(bool is_mc,TString era);

  double GetEETriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta);
  double GetMMTriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta);
  double GetEMTriggerScaleFactor(Float_t l1_pt, Float_t l2_pt, Float_t l1_eta, Float_t l2_eta);
  NanoEffCorrection_t GetChargeFlipScaleFactor(TString File_name, TString Fit_type, Float_t pt1, Float_t eta1, Int_t charge1, Float_t pt2, Float_t eta2, Int_t charge2);

  ~NanoScaleFactorsWrapper();

 private:
  std::map<TString,TH2 *> scaleFactorsH_;
  int era_;
  Float_t chargeflip_SF,chargeflip_evwgt_woSF,chargeflip_evwgt_SF;
  int chargeflip_count, chargeflip_ratio_count;

};

#endif
