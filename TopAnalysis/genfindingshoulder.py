import argparse
import ROOT
import json
import numpy as np

def lepton_order(t):
  if(t.nl<2): return[-1,-1]
  l1_index = -1
  l2_index = -1
  l1_pt = -1.
  l2_pt = -1.
  for i in range(t.nl):
    if(t.l_pt[i]>l1_pt):
      l2_index = l1_index
      l2_pt = l1_pt
      l1_index = i
      l1_pt = t.l_pt[i]
    elif (t.l_pt[i]>l2_pt):
      l2_index = i
      l2_pt = t.l_pt[i]
  return [l1_index,l2_index]
      
def invM(t):
  if(t.nl<2): return 0
  e = lepton_order(t)
  pt1 = t.l_pt[e[0]]
  pt2 = t.l_pt[e[1]]
  eta1 = t.l_eta[e[0]]
  eta2 = t.l_eta[e[1]]
  phi1 = t.l_phi[e[0]]
  phi2 = t.l_phi[e[1]]
  return round((2*pt1*pt2*(np.cosh(eta1-eta2)-np.cos(phi1-phi2)))**0.5,1)

def geninvM(t):
  if(t.nl<2): return 0
  e = lepton_order(t)
  g = [t.l_g[e[0]],t.l_g[e[1]]]
  if(g[0]<0 or g[1]<0): return 0
  pt1 = t.g_pt[g[0]]
  pt2 = t.g_pt[g[1]]
  eta1 = t.g_eta[g[0]]
  eta2 = t.g_eta[g[1]]
  phi1 = t.g_phi[g[0]]
  phi2 = t.g_phi[g[1]]
  return round((2*pt1*pt2*(np.cosh(eta1-eta2)-np.cos(phi1-phi2)))**0.5,1)


#Setting Parse
parser = argparse.ArgumentParser(description="Process the root file")
parser.add_argument('fin', type=str)
args = parser.parse_args()

#Read File
inf = ROOT.TFile.Open(args.fin)
t = inf.Get("analysis/data")
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
    Zmass = invM(t)
    genZmass = geninvM(t)
    if ( lb < Zmass and Zmass < ub ):
      counting+=1
      e = lepton_order(t)
      g = [t.l_g[e[0]],t.l_g[e[1]]]
      print("Entry: "+str(i))
      print("--[l1]--  pt: "+str(round(t.l_pt[e[0]],1))+" eta: "+str(round(t.l_eta[e[0]],2))+" phi: "+str(round(t.l_phi[e[0]],2))+" Zmass: "+str(Zmass))
      if(g[0]<0): print("--[g1]-- not matching")
      else: print("--[g1]--  pt: "+str(round(t.g_pt[g[0]],1))+" eta: "+str(round(t.g_eta[g[0]],2))+" phi: "+str(round(t.g_phi[e[0]],2))+" Zmass: "+str(genZmass))
      print("--[l2]--  pt: "+str(round(t.l_pt[e[1]],1))+" eta: "+str(round(t.l_eta[e[1]],2))+" phi: "+str(round(t.l_phi[e[1]],2))+" Zmass: "+str(Zmass))
      if(g[1]<0): print("--[g2]-- not matching")
      else: print("--[g2]--  pt: "+str(round(t.g_pt[g[1]],1))+" eta: "+str(round(t.g_eta[g[1]],2))+" phi: "+str(round(t.g_phi[g[1]],2))+" Zmass: "+str(genZmass))

