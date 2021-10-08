import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np

plotdir = "/eos/user/t/tihsu/DataAnalysis/Charge_Flip_woSB_genStudy/"
apply_SF = False

# Read Files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
AllFiles = os.listdir(targetdir)

if(apply_SF):
  sigdir = os.path.join(workdir+"/../analysis_2017_chargeflip_SF_combine/Chunks")
  sigfile = os.listdir(sigdir)
# Read json Files
jsonfile = open(os.path.join(workdir+"/../samples_2017.json"))
samplesList = json.load(jsonfile, encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()



# Histograms
h_oc_Zmass_MC = ROOT.TH1F('oc_Zmass_MC',"; M(l,l) [GeV] ; Events",150,0,150)
h_oc_Zmass_data = ROOT.TH1F('oc_Zmass_data',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_MC = ROOT.TH1F('ss_Zmass_MC',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_data = ROOT.TH1F('ss_Zmass_data',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_DMinC = ROOT.TH1F('ss_Zmass_DMinC',"; M(l,l) [GeV] ; Events",150,0,150)

diagrams = ['ee_ss_Zmass', 'ee_oc_Zmass']
Zmass_result = dict()
sig_list = ["DY","t#bar{t}"]


for s,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in AllFiles:
    if s in ff:
      inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
      try:
        h_oc = inf.Get('ee_oc_Zmass')
        h_ss = inf.Get('ee_ss_Zmass')
        if ("MC" in ff and (desc[3] in sig_list) and not apply_SF):
          print("signal"+ff)
          h_oc_Zmass_MC.Add(h_oc_Zmass_MC,h_oc,1.0,weight)
          h_ss_Zmass_MC.Add(h_ss_Zmass_MC,h_ss,1.0,weight)
        elif("Data" in ff):
          h_oc_Zmass_data.Add(h_oc_Zmass_data,h_oc,1.0,1.0)
          h_ss_Zmass_data.Add(h_ss_Zmass_data,h_ss,1.0,1.0)
          #else:
          #  h_oc_Zmass_data.Add(h_oc_Zmass_data,h_oc,1.0,-1.*weight)
          #  h_ss_Zmass_data.Add(h_ss_Zmass_data,h_ss,1.0,-1.*weight)
          print(ff)
      except:
        print(ff+":no histogram")
        pass
      inf.Close()

print("-------------------------")
if(apply_SF):
  for s,desc in samplesList:
    weight = desc[0]*41.5*1000
    for ff in sigfile:
      if s in ff:
        inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
        try:
          h_oc = inf.Get('ee_oc_Zmass')
          h_ss = inf.Get('ee_ss_Zmass')
          if ("MC" in ff and (desc[3] in sig_list)):
            print("signal")
            h_oc_Zmass_MC.Add(h_oc_Zmass_MC,h_oc,1.0,weight)
            h_ss_Zmass_MC.Add(h_ss_Zmass_MC,h_ss,1.0,weight)
        except:
          print("no histogram")
          pass
        print(ff)
        inf.Close()

c = ROOT.TCanvas('c','c',0,0,600,400)
ROOT.gStyle.SetOptStat(0000)

numberlist = [0.,0.,0.,0.]
numberlist[0] = h_ss_Zmass_data.Integral()
numberlist[1] = h_ss_Zmass_MC.Integral()
numberlist[2] = h_oc_Zmass_data.Integral()
numberlist[3] = h_oc_Zmass_MC.Integral()

MC_SF = (numberlist[0]+numberlist[2])/(numberlist[1]+numberlist[3])


h_oc_Zmass_MC.GetXaxis().SetRange(70,110)
h_oc_Zmass_MC.Scale(MC_SF)
h_oc_Zmass_MC.Draw("HIST")
c.Update()

h_oc_Zmass_data.GetXaxis().SetRange(70,110)
h_oc_Zmass_data.SetMarkerStyle(20)
h_oc_Zmass_data.Draw("E2 SAME")
c.Update()

print(h_ss_Zmass_MC.GetMaximum())

rightmax = 1.1 * float(h_ss_Zmass_MC.GetMaximum())
scale = ROOT.gPad.GetUymax()/rightmax

h_ss_Zmass_MC.GetXaxis().SetRange(70,110)
h_ss_Zmass_MC.SetLineColor(2)
h_ss_Zmass_MC.Scale(scale*MC_SF)
h_ss_Zmass_MC.Draw("HIST SAME")
c.Update()

h_ss_Zmass_data.GetXaxis().SetRange(70,110)
h_ss_Zmass_data.SetLineColor(2)
h_ss_Zmass_data.Scale(scale)
h_ss_Zmass_data.SetMarkerStyle(20)
h_ss_Zmass_data.SetMarkerColor(2)
h_ss_Zmass_data.Draw("E2 SAME")
c.Update()
legend = ROOT.TLegend(0.1,0.7,0.44448,0.9)
legend.AddEntry(h_ss_Zmass_data,"Data(SS)","p")
legend.AddEntry(h_ss_Zmass_MC,"MC(SS)","l")
legend.AddEntry(h_oc_Zmass_data,"Data(OC)","p")
legend.AddEntry(h_oc_Zmass_MC,"MC(OC)","l")
legend.Draw("SAME")

axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(),ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymax(),0,rightmax,20,"+L")
axis.SetLineColor(2)
axis.SetLabelColor(2)
axis.Draw()
c.Update()

c.SaveAs(plotdir+'Zmass.png')
c.SaveAs(plotdir+'Zmass.pdf')


# Fitting
w = ROOT.RooWorkspace()
#w.factory('Gaussian::pdf(x[0,150],mean[90,50,130],gamma[0,20])')
w.factory('BreitWigner::BW(x[50,130],m[90,70,110],#Gamma[0,8.])')
w.factory('Exponential::e(x,alpha[-5,-10,10])')
w.factory('SUM::pdf(a[9000,0,100000000]*BW,b[10,0,100000000]*e)')
w.factory('Gaussian::g(x,M[90,70,110],sigma[0,20])')
#w.factory('CBShape::CB(x,m[0,150],s[0,150],a[0,150],n[-5,5])')
#w.factory('FCONV::pdf(x,BW,CB)')

pdf = w.pdf('pdf')
x = w.var('x')
gamma = w.var('#Gamma')
mean = w.var('m')
BW = w.pdf('BW')
g = w.pdf('g')


ROOT.gStyle.SetOptStat(1101)
ROOT.gStyle.SetOptFit(1011);

meanlist = [0,0,0,0]
gammalist = [0,0,0,0]
widthlist = [0,0,0,0]
wglist = [0,0,0,0]
goallist = [0,0,0,0]
sig_ratio = 0.935
name_list = ["SS_data","SS_MC","OS_data","OS_MC"]
h_list = [h_ss_Zmass_data,h_ss_Zmass_MC,h_oc_Zmass_data,h_oc_Zmass_data]

for i in range(4):
  dh_ss_data = ROOT.RooDataHist('dh_'+name_list[i],'dh_'+name_list[i],ROOT.RooArgList(x),h_list[i])
  pdf.fitTo(dh_ss_data)
  xframe = x.frame()
  xframe.SetTitle(name_list[i])
  xframe.GetXaxis().SetTitle("M(l,l) [GeV]")
  xframe.GetYaxis().SetTitle("Events")
  dh_ss_data.plotOn(xframe)
  pdf.plotOn(xframe)
#  xframe.Draw()
  pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
#  xframe.addObject(pt)
  #xframe.Draw()
  BW.paramOn(xframe, ROOT.RooFit.Layout(0.55));
  xframe.Draw()
  c.Update()
  mean.Print()
  gamma.Print()
  c.SaveAs(plotdir+'fit_'+name_list[i]+'.png')
  c.SaveAs(plotdir+'fit_'+name_list[i]+'.pdf')

  meanlist[i] = mean.getValV()
  gammalist[i] = gamma.getValV()
  ig = BW.createIntegral(ROOT.RooArgSet(x))
  goal = 0.
  width = 0.
  while(goal<sig_ratio):
    width+=0.0001
    x.setRange(mean.getValV()-width,mean.getValV()+width)
    goal = ig.getVal()
    x.setRange(50,130)
    goal = goal/ig.getVal()
  print(str(width)+" : "+str(goal)+" "+ str(width/gamma.getValV()))
  widthlist[i] = width
  goallist[i] = goal
  wglist[i] = width/gamma.getValV()
#numberlist[0] = h_ss_Zmass_data.Integral()

#numberlist[3] = h_oc_Zmass_MC.Integral()
#ig.PlotOn(xframe)
#c.Update()
#c.SaveAs(plotdir+"CDF.png")
meanlist[0] = meanlist[1]
meanlist[2] = meanlist[3]
widthlist[0] = widthlist[1]
widthlist[2] = widthlist[3]
print("ss data ; ss MC ; oc data; oc MC")
print("mean ; gamma ; width ; goal ; number")
print(meanlist)
print(gammalist)
print(widthlist)
print(goallist)
print(wglist)
print(numberlist)
