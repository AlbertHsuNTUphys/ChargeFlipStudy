import argparse
import ROOT
import json
import numpy as np

def invM(t):
  pt1 = t.t_pt_l1
  pt2 = t.t_pt_l2
  eta1 = t.t_eta_l1
  eta2 = t.t_eta_l2
  phi1 = t.t_phi_l1
  phi2 = t.t_phi_l2
  return round((2*pt1*pt2*(np.cosh(eta1-eta2)-np.cos(phi1-phi2)))**0.5,1)

#Setting Parse
parser = argparse.ArgumentParser(description="Process the root file")
parser.add_argument('fin', type=str)
args = parser.parse_args()

#Read File
inf = ROOT.TFile.Open(args.fin)
t = inf.Get("TreeInput")
print(t)

#Main Code
while(1):
  number = input('number of events :')
  lb = input('Lower Bound :')
  ub = input('Upper Bound :')
  counting = 0
  for i in range(t.GetEntriesFast()):
    t.GetEntry(i)
    if (counting == number): continue
    if ( lb < t.t_Zmass and t.t_Zmass < ub ):
      counting+=1
      print("Entry: "+str(i))
      print("--[l1]--  pt: "+str(round(t.t_pt_l1,1))+" eta: "+str(round(t.t_eta_l1,2))+" phi: "+str(round(t.t_phi_l1,2))+" Zmass: "+str(round(t.t_Zmass,1))+"/"+str(invM(t)))
      print("--[l2]--  pt: "+str(round(t.t_pt_l2,1))+" eta: "+str(round(t.t_eta_l2,2))+" phi: "+str(round(t.t_phi_l2,2))+" Zmass: "+str(round(t.t_Zmass,1))+"/"+str(invM(t)))
