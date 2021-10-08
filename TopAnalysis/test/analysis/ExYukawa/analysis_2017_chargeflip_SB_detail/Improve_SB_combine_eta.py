import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import math

plotdir = "/eos/user/t/tihsu/DataAnalysis/Charge_Flip_woSB/"
tag = 'combine'
ori_detail = True
# Read Files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)

# Read json Files
jsonfile = open(os.path.join(workdir+"/../samples_2017.json"))
samplesList = json.load(jsonfile,encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Basic Setting
pt_region = np.array(( 20.0, 50.0, 100., 200., 300.))
eta_region = np.array((0.0, 0.8, 1.479, 2.4))
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
channel = ['MC_oc','MC_ss','data_oc','data_ss','back_oc','back_ss']
Number = dict()
Pij = dict()
Vij = dict()
P_fit = dict()

N_MC = 0.
N_data = 0.

for ch in channel:
  Number[ch] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]

skip_list = ["MC13TeV_2017_QCDEM_120to170_2.root", "MC13TeV_2017_QCDEM_120to170_1.root", "MC13TeV_2017_QCDEM_80to120_1.root","MC13TeV_2017_QCDEM_80to120_3.root","MC13TeV_2017_QCDEM_120to170_4.root","MC13TeV_2017_WJets_mlm_49.root"]
# Run Through All files
for s ,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in Allfile:
    if s in ff:
      try:
        inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
        try:
          t = inf.Get("TreeInput")
        except:
          t = inf.Get("Chunks/"+ff+"/TreeInput")
        n_entry = t.GetEntries()
        if "MC" in s and (desc[3] == "DY" or desc[3]=="t#bar{t}"):
          print("Signal: "+ff)
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  Number['MC_ss'][i][j][ii][jj] += t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
                  Number['MC_oc'][i][j][ii][jj] += t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
                  if (math.isnan(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):
                    print("******* "+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+str(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]))
                  N_MC += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight + float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight
        elif "Data" in s:
          print("Data: "+ff)
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  Number['data_ss'][i][j][ii][jj] += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  Number['data_oc'][i][j][ii][jj] += float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  N_data += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]) + float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  if (math.isnan(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):
                    print("******* "+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+str(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]))
        elif(ff not in skip_list):
          print("Background: "+ff)
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  #if (math.isnan(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):continue
                  #if (math.isnan(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):continue
                  Number['back_ss'][i][j][ii][jj] += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight
                  Number['back_oc'][i][j][ii][jj] += float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight
                  if (math.isnan(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):
                    print("******* "+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+str(float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])))
                  N_MC += float(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight + float(t.t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])*weight
        inf.Close()
      except:
        print(ff+"-->Trigger exception")
print(N_data)
print(N_MC)
normscale = 1.
"""for i in range(pt_bins):
  for j in range(eta_bins):
    for ii in range(pt_bins):
      for jj in range(eta_bins):
        N_MC += Number['MC_ss'][i][j][ii][jj]
        N_MC += Number['MC_oc'][i][j][ii][jj]
        N_data += Number['data_ss'][i][j][ii][jj]
        N_data += Number['data_oc'][i][j][ii][jj]
        N_MC += Number['back_ss'][i][j][ii][jj]
        N_MC += Number['back_oc'][i][j][ii][jj]
"""
print(N_data)
print(N_MC)
normscale = N_data/N_MC
print(normscale)
normscale = 1.
for i in range(pt_bins):
  for j in range(eta_bins):
    for ii in range(pt_bins):
      for jj in range(eta_bins):
        Number['MC_ss'][i][j][ii][jj]*=normscale
        Number['MC_oc'][i][j][ii][jj]*=normscale
        Number['data_ss'][i][j][ii][jj]-=Number['back_ss'][i][j][ii][jj]*normscale
        Number['data_oc'][i][j][ii][jj]-=Number['back_oc'][i][j][ii][jj]*normscale

# Combine symmetric result
for h in Number:
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          if(ii>i or jj>j):
            Number[h][ii][jj][i][j]+=Number[h][i][j][ii][jj]
            Number[h][i][j][ii][jj] = 0

# Combine 100-200GeV with 200 GeV up
if(tag=='combine' and ori_detail):
  pt_region = np.array(( 20.0, 50.0, 100.0,200., 300.))
  eta_region = np.array((0.0, 1.479, 2.4))
  pt_bins = len(pt_region)-1
  eta_bins = len(eta_region)-1
  for h in Number:
   for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          if(j==eta_bins-1 and jj==eta_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i][j+1][ii][jj+1]
            Number[h][i][j+1][ii][jj+1] = 0
          if(j==eta_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i][j+1][ii][jj]
            Number[h][i][j][ii][jj+1] = 0
          if(ii==pt_bins-1):
            Number[h][i][j][ii][jj] += Number[h][i][j][ii][jj+1]
            Number[h][i][j][ii+1][jj] = 0

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
for i in range(pt_bins):
  for j in range(eta_bins):
    for ii in range(pt_bins):
      for jj in range(eta_bins):
        print(str(i)+str(j)+str(ii)+str(jj)+" : data: "+str(Pij['data'][i][j][ii][jj])+" MC: "+str(Pij['MC'][i][j][ii][jj]))
        print(str(Number['data_ss'][i][j][ii][jj])+" "+str(Number['MC_ss'][i][j][ii][jj]))


# Use chi2 to fit
    
for h in DvMC:
  print(Pij[h])
  gMinuit = ROOT.TMinuit(pt_bins*eta_bins)

  def fcn(npar, gin, f, par, iflag):
      chi2 = 0.  
      Likelihood = 0.
      for i in range(pt_bins):
        for j in range(eta_bins):
          for ii in range(pt_bins):
            for jj in range(eta_bins):
              if (not Vij[h][i][j][ii][jj]==0.) and Number[h+'_oc'][i][j][ii][jj]>30 and Number[h+'_ss'][i][j][ii][jj]>1:
                P1 = par[i*eta_bins+j]
                P2 = par[ii*eta_bins+jj]
                N = int(round(Number[h+'_oc'][i][j][ii][jj]+Number[h+'_ss'][i][j][ii][jj]))
                Nsc = int(round(Number[h+'_ss'][i][j][ii][jj]))
                chi2 += ((Pij[h][i][j][ii][jj]-(P1+P2-P1*P2))**2)/(Vij[h][i][j][ii][jj])
                Likelihood += Nsc*np.math.log((N*(P1+P2)))-N*(P1+P2)-np.math.log(np.math.factorial(Nsc))

      f[0] = chi2

  gMinuit.SetFCN(fcn)
  for i in range(pt_bins):
    for j in range(eta_bins):
      init_val = 0.0001
      if Number[h+'_oc'][i][j][i][j]>5000:
        init_val = 1.-(1.-Pij[h][i][j][i][j])**0.5
        print(init_val)
      gMinuit.DefineParameter(i*eta_bins+j,"P"+str(i)+str(j),init_val,0.00000001,0.,0.1) 
  gMinuit.Command("Minuit2")
#  r = gMinuit.save()
# Result Plot
  c = ROOT.TCanvas(h+'c','c',600,600)
  ROOT.gStyle.SetOptStat("kFALSE")
  ROOT.gStyle.SetPaintTextFormat(".2e");
  ROOT.gStyle.SetPalette(69);
#  r.correlationHist(h).Draw('colz')
  c.Update()
  c.SaveAs(plotdir+h+'_CorrelationHist_SB_combine_Cov.png')
  
  h_chargeflip = ROOT.TH2D(h+"CFRate",";P_{T}[GeV] ; \eta;"+h+"_ChargeFlip rate",pt_bins,pt_region,eta_bins,eta_region)

  result = [[0. for j in range(eta_bins)] for i in range(pt_bins)]
  error = [[0. for j in range(eta_bins)] for i in range(pt_bins)]

  for i in range(pt_bins):
    for j in range(eta_bins):
      result = ROOT.double(0.)
      error = ROOT.double(0.)
      gMinuit.GetParameter(i*eta_bins+j,result,error)
      P_fit[h][i][j] = result
      h_chargeflip.SetBinContent(i+1,j+1,result)
      h_chargeflip.SetBinError(i+1,j+1,error)
  if h=='data':
    h_data = h_chargeflip.Clone()
  else:
    h_MC = h_chargeflip.Clone()
  c.SetLogx()
  h_chargeflip.Draw('COLZTEXT e')
  c.Update()
  c.SaveAs(plotdir+h+'_CFRate_SB_'+tag+'_Cov.png')
h_data.Divide(h_data,h_MC)
h_data.Draw('COLZTEXT e')
c.Update()
c.SaveAs(plotdir+'CFRate_MDRatio_SB_'+tag+'_Cov.png')

for i in range(pt_bins):
  for j in range(eta_bins):
    SF[i][j] = P_fit['data'][i][j]/P_fit['MC'][i][j]

print(P_fit['MC'])
print(SF)
