import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import math

plotdir = "/eos/user/t/tihsu/DataAnalysis/Charge_Flip_woSB_genStudy/"


# Read Files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)

# Read json Files
jsonfile = open(os.path.join(workdir+"/../samples_2017_chargeflip.json"))
samplesList = json.load(jsonfile,encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Basic Setting
pt_region = np.array(( 20.0, 50.0, 100.0, 300.))
eta_region = np.array((0.0, 0.8, 1.479, 2.4))
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
channel = ['gen_oc','gen_ss','gen_mu_oc','gen_mu_ss']
Number = dict()
Pij = dict()
Vij = dict()

for ch in channel:
  Number[ch] = [[0. for i in range(eta_bins)] for j in range(pt_bins)]
#Mean = [90.620,89.991,90.218,89.113]
#Gamma = [4.267,5.230,4.691,6.350]
#Zmass = dict()
#Zwidth = dict()
#for i in range(len(channel)):
#  Zmass[channel[i]] = Mean[i]
#  Zwidth[channel[i]] = Gamma[i]

# Run Through All files
for s ,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in Allfile:
    if s in ff:
      print(ff)
      try:
        inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
        t = inf.Get("TreeInput")
        n_entry = t.GetEntries()
        if "MC" in s:
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
                  Number['gen_ss'][i][j] += t.t_g_Nss[j + eta_bins*i]*weight
                  Number['gen_oc'][i][j] += t.t_g_Noc[j + eta_bins*i]*weight
                  Number['gen_mu_oc'][i][j] += t.t_gm_Noc[j + eta_bins*i]*weight
                  Number['gen_mu_ss'][i][j] += t.t_gm_Nss[j + eta_bins*i]*weight
        inf.Close()
      except:
        print(ff+"-->Trigger exception")
# Combine 100-200GeV with 200GeV UP
pt_region = np.array(( 20.0, 50.0, 100.0, 300.))
eta_region = np.array((0.0, 0.8, 1.479, 2.4))
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
"""for h in Number:
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          if(i==pt_bins-1 and ii==pt_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i+1][j][ii+1][jj]
            Number[h][i+1][j][ii+1][jj] = 0
          if(i==pt_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i+1][j][ii][jj]
            Number[h][i+1][j][ii][jj] = 0
          if(ii==pt_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i][j][ii+1][jj]
            Number[h][i][j][ii+1][jj] = 0
          if(ii>i or jj>j):
            Number[h][ii][jj][i][j]+=Number[h][i][j][ii][jj]
            Number[h][i][j][ii][jj] = 0"""

# Calculate Probability
DvMC = ['gen', 'gen_mu']
for h in DvMC:
  Pij[h] = [[0. for i in range(eta_bins)] for j in range(pt_bins)]
  Vij[h] = [[0. for i in range(eta_bins)] for j in range(pt_bins)]
  for i in range(pt_bins):
    for j in range(eta_bins):
#          print("l1_pt = " + str(pt_region[i]) + " l1_eta = " + str(eta_region[j]) + " l2_pt = " + str(pt_region[ii]) + " l2_eta = " + str(eta_region[jj]))
#          print('ss : '+str(Number[h+'_ss'][i][j][ii][jj])+' oc : '+str(Number[h+'_oc'][i][j][ii][jj]))
      if(1):
        N_ss = Number[h+'_ss'][i][j]
        N_T = N_ss + Number[h+'_oc'][i][j]
        if(N_T>0. and N_ss>0.):
          Pij[h][i][j] = N_ss/N_T 
          Vij[h][i][j] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3)-2*(N_ss**1.5)/(N_T**2.5)
#            Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3) # without consider covariance 
        else:
          print("raise Number counting error")
          pass


# Use chi2 to fit
Fout = ROOT.TFile.Open("genCFrate.root","RECREATE")    
for h in DvMC:
  print(Pij[h])
#  r = gMinuit.save()
# Result Plot
  c = ROOT.TCanvas(h+'c','c',600,600)
  ROOT.gStyle.SetOptStat("kFALSE")
  ROOT.gStyle.SetPaintTextFormat(".2e");
  ROOT.gStyle.SetPalette(69);
  c.SetRightMargin(0.15);
  c.SetTopMargin(0.15);

#  r.correlationHist(h).Draw('colz')
  c.Update()
#  c.SaveAs(plotdir+h+'_CorrelationHist_SB_combine_Cov.png')
  
  h_chargeflip = ROOT.TH2D(h+"CFRate",";P_{T}[GeV] ; |\eta|;",pt_bins,pt_region,eta_bins,eta_region)

#  result = [[0. for j in range(eta_bins)] for i in range(pt_bins)]
#  error = [[0. for j in range(eta_bins)] for i in range(pt_bins)]

  for i in range(pt_bins):
    for j in range(eta_bins):
      h_chargeflip.SetBinContent(i+1,j+1,Pij[h][i][j])
      h_chargeflip.SetBinError(i+1,j+1,Vij[h][i][j]**0.5)
  c.SetLogx()
  c.SetLogz()
  h_chargeflip.Draw('COLZTEXT e')
  c.Update()
  c.SaveAs(plotdir+h+'_CFRate_woSB_combine_Cov.png')
  c.SaveAs(plotdir+h+'_CFRate_woSB_combine_Cov.pdf')
  h_chargeflip.Write()

Fout.Close()
