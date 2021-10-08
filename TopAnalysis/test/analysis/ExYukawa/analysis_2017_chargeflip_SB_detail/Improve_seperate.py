import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import math
from tempfile import TemporaryFile

plotdir = "/eos/user/t/tihsu/DataAnalysis/Charge_Flip_woSB/"


# Read Files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)

# Read json Files
jsonfile = open(os.path.join(workdir+"/../samples_2017_chargeflip.json"))
samplesList = json.load(jsonfile,encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Basic Setting
pt_region = np.array((20.,30.,40., 50.,60.,70.,80.,90., 100., 120.,150.,200., 300.));
eta_region = np.array((0., 0.5, 0.7, 0.75, 1.0, 1.2, 1.479, 1.7, 1.9, 2.1, 2.3, 2.4));
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
channel = ['MC_oc','MC_ss','data_oc','data_ss']
Number = dict()
Pij = dict()
Vij = dict()
P_fit = dict()

for ch in channel:
  Number[ch] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
Mean = [90.620,89.991,90.218,89.113]
Gamma = [4.267,5.230,4.691,6.350]
Zmass = dict()
Zwidth = dict()
for i in range(len(channel)):
  Zmass[channel[i]] = Mean[i]
  Zwidth[channel[i]] = Gamma[i]

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
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  Number['MC_ss'][i][j][ii][jj] += t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
                  Number['MC_oc'][i][j][ii][jj] += t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
        if "Data" in s:
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  Number['data_ss'][i][j][ii][jj] += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  Number['data_oc'][i][j][ii][jj] += float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  if (math.isnan(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):
                    print("******* "+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+str(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]))
        inf.Close()
      except:
        print(ff+"-->Trigger exception")
# Combine 100-200GeV with 200GeV UP
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

P_fit['data'] = [[0. for i in range(eta_bins)] for j in range(pt_bins)]
P_fit['MC'] = [[0. for i in range(eta_bins)] for j in range(pt_bins)]
SF = [[0. for i in range(eta_bins)] for j in range(pt_bins)]

# Calculate Probability
DvMC = ['data','MC']
for h in DvMC:
  Pij[h] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
  Vij[h] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
#          print("l1_pt = " + str(pt_region[i]) + " l1_eta = " + str(eta_region[j]) + " l2_pt = " + str(pt_region[ii]) + " l2_eta = " + str(eta_region[jj]))
#          print('ss : '+str(Number[h+'_ss'][i][j][ii][jj])+' oc : '+str(Number[h+'_oc'][i][j][ii][jj]))
          if(1):
            N_ss = Number[h+'_ss'][i][j][ii][jj]
            N_T = N_ss + Number[h+'_oc'][i][j][ii][jj]     
            if(N_T>0. and N_ss>0.):
              Pij[h][i][j][ii][jj] = N_ss/N_T
              Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3)-2*(N_ss**1.5)/(N_T**2.5)
#            Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3) # without consider covariance 
          else:
            print("raise Number counting error")
            pass


with open('seperate_CFrate.npy', 'wb') as f:
    np.save(f, np.array(Pij['data']))
    np.save(f, np.array(Pij['MC']))
    np.save(f,np.array(Vij['data']))
    np.save(f,np.array(Vij['MC']))

# Use chi2 to fit
c = ROOT.TCanvas('c','c',600,600)
ROOT.gStyle.SetOptStat('kFALSE')
ROOT.gStyle.SetMarkerStyle(20)
c.Update()
    
for h in DvMC:

  print(Pij[h])
  gMinuit = ROOT.TMinuit(pt_bins+eta_bins)

  def fcn(npar, gin, f, par, iflag):
      chi2 = 0.  
      for i in range(pt_bins):
        for j in range(eta_bins):
          for ii in range(pt_bins):
            for jj in range(eta_bins):
              if (not Vij[h][i][j][ii][jj]==0.) and Number[h+'_oc'][i][j][ii][jj]>0 and Number[h+'_ss'][i][j][ii][jj]>0:
                P1 = par[i+eta_bins]*par[j]
                P2 = par[ii+eta_bins]*par[jj]
                chi2 += (Pij[h][i][j][ii][jj]-(P1+P2-P1*P2))**2/(Vij[h][i][j][ii][jj])
      f[0] =  chi2 

  gMinuit.SetFCN(fcn)
  for i in range(pt_bins):
    gMinuit.DefineParameter(eta_bins+i,"P"+str(i),0.01,0.0001,0.0,10.0)
  for j in range(eta_bins):
    gMinuit.DefineParameter(j,"eta"+str(j),0.01,0.0001,0.0,1.0)
  gMinuit.Command("Minuit2")
  
  #h_chargeflip = ROOT.TH2D(h+"CFRate",";P_{T}[GeV] ; \eta;"+h+"_ChargeFlip rate",pt_bins,pt_region,eta_bins,eta_region)
  h_E_plot = ROOT.TH1F(h+"CF_E",";|\eta| ; CF Rate;"+h+"_ChargeFlip rate",eta_bins,eta_region)
  h_P_plot = ROOT.TH1F(h+"CF_P",";PT[GeV] ; CF Rate;"+h+"_ChargeFlip rate",pt_bins,pt_region)

  for i in range(pt_bins):
    result = ROOT.double(0.)
    error = ROOT.double(0.)
    gMinuit.GetParameter(eta_bins+i,result,error)
    h_P_plot.SetBinContent(i+1,result)
    h_P_plot.SetBinError(i+1,error)
  for i in range(eta_bins):
    result = ROOT.double(0.)
    error = ROOT.double(0.)
    gMinuit.GetParameter(i,result,error)
    h_E_plot.SetBinContent(i+1,result)
    h_E_plot.SetBinError(i+1,error)

  if h=='data':
    h_data_E = h_E_plot.Clone()
    h_data_P = h_P_plot.Clone()
  else:
    h_MC_E = h_E_plot.Clone()
    h_MC_P = h_P_plot.Clone()

c.SetLogx()
c.SetLogy()
h_data_E.GetXaxis().SetRangeUser(0,2.4)
h_MC_E.GetXaxis().SetRangeUser(0,2.4)
h_data_E.Draw("E2")
h_MC_E.SetMarkerColor(2)
h_MC_E.Draw("E2 SAME")
c.SaveAs(plotdir+'P(Eta)_chargeflip.png')
h_data_P.Draw("E2")
h_MC_P.SetMarkerColor(2)
h_MC_P.Draw("E2 SAME")
c.SaveAs(plotdir+'P(Pt)_chargeflip.png')

rp_E = ROOT.TRatioPlot(h_data_E, h_MC_E);
rp_E.GetLowerRefYaxis().SetTitle("Data/MC")
c.SetTicks(0, 1);
rp_E.Draw();
rp_E.GetLowYaxis().SetNdivisions(505);
c.SaveAs(plotdir+"P(Eta)_ratio.png")

rp_P = ROOT.TRatioPlot(h_data_P, h_MC_P);
rp_P.GetLowerRefYaxis().SetTitle("Data/MC")
c.SetTicks(0, 1);
rp_P.Draw();
rp_P.GetLowYaxis().SetNdivisions(505);
c.SaveAs(plotdir+"P(Pt)_ratio.png")





