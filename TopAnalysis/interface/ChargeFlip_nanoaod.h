#ifndef _chargeflip_nanoaod_h_
#define _chargeflip_nanoaod_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/NanoSelectionTool.h"

void RunChargeFlip_nanoaod(const TString in_fname,
                                 TString outname,
                                 TH1F   *normH,
                                 TString era,
                                 Bool_t  debug);

#endif
