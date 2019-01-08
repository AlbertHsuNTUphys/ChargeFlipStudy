#!/usr/bin/env python
import os,sys
import ROOT
from optparse import OptionParser
from collections import defaultdict
import numpy as np

COLORS=[1, ROOT.kOrange,  ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]

def doNuisanceReport(args,outdir):

    """compare postfit nuisances"""

    blackList=['seed','itoy','status','errstatus','scanidx','edmval','nllval','nllvalfull','dnllval','chisq','chisqpartial','ndof','ndofpartial']

    allVars=[]
    results=[]
    for i in xrange(0,len(args)):

        title,url=args[i].split(':')

        inF=ROOT.TFile.Open(url)
        fitres=inF.Get('fitresults')

        #get the names of systematics
        varVals=defaultdict(dict)
        fitres.GetEntry(0)
        for b in fitres.GetListOfBranches():
            bname=b.GetName()
            if bname in blackList: continue
            if 'ttbar_mu' in bname: continue
            if '_In' in bname : continue
            if '_gen' in bname : continue
            if '_minos' in bname : continue
            if '_err' in bname : continue
            varVals[bname]=[]
            
        ntoys=fitres.GetEntriesFast()

        #if several toys compute the median and the 68% CL
        if ntoys>1:
            for n in range(ntoys):
                fitres.GetEntry(n)
                for v in varVals: 
                    varVals[v].append( getattr(fitres,v) )
            
            for v in varVals:
                qtl=np.percentile(varVals[v],[16,50,84])
                varVals[v]=[qtl[1],qtl[1]-qtl[0],qtl[2]-qtl[1]]

        #use post-fit uncertainty for asimov or data fit 
        else:

            fitres.GetEntry(0)
            for v in varVals: 
                varVals[v]=[ getattr(fitres,v),getattr(fitres,v+'_err'),getattr(fitres,v+'_err') ]

        for v in varVals: allVars.append(v)
        results.append( (title,varVals) )

    allVars=list(set(allVars))
    

    #show nuisances
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.3)
    c.SetTopMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(1.0)
    c.SetGridy(True)
    npergroup=20
    ngroups=len(allVars)/npergroup+1
    for ig in range(ngroups):
        first=npergroup*ig
        last=min(npergroup*(ig+1),len(allVars))
        varList=allVars[first:last]

        #prepare a frame
        npars=len(varList)
        frame=ROOT.TH2F('frame',';#hat{#theta};Nuisance parameter',1,-3,3,npars,0,npars)
        frame.SetDirectory(0)
        for ybin in range(npars):
            frame.GetYaxis().SetBinLabel(ybin+1,'#color[%d]{%s}'%((ybin%2)*10+1,varList[ybin]))
        frame.Draw()
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(0.8)
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().SetTitleOffset(3.2)

        #1-sigma
        gr1s=ROOT.TGraph()
        gr1s.SetName('gr1s')
        gr1s.SetMarkerStyle(1)
        gr1s.SetMarkerColor(19)
        gr1s.SetLineColor(19)
        gr1s.SetFillStyle(1001)
        gr1s.SetFillColor(19)
        gr1s.SetPoint(0,-1,0)
        gr1s.SetPoint(1,-1,npars)
        gr1s.SetPoint(2,1,npars)
        gr1s.SetPoint(3,1,0)
        gr1s.SetPoint(4,-1,0)

        #2-sigma
        gr2s=gr1s.Clone('gr2s')
        gr2s.SetMarkerColor(18)
        gr2s.SetLineColor(18)
        gr2s.SetFillStyle(1001)
        gr2s.SetFillColor(18)
        gr2s.SetPoint(0,-2,0)
        gr2s.SetPoint(1,-2,npars)
        gr2s.SetPoint(2,2,npars)
        gr2s.SetPoint(3,2,0)
        gr2s.SetPoint(4,-2,0)

        gr2s.Draw('f')
        gr1s.Draw('f')

        txt=ROOT.TLatex()
        txt.SetTextFont(42)
        txt.SetTextSize(0.03)
        txt.SetTextColor(ROOT.kGray+3)
        for delta,title in [(1.0,'-1#sigma'),(2,'+2#sigma'),(-1,'-1#sigma'),(-2,'-2#sigma')]:
            txt.DrawLatex(delta-0.2,frame.GetYaxis().GetXmax()+0.2,title)

        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.05,0.955,'#bf{CMS} #it{Preliminary}')
        txt.SetTextAlign(31)
        txt.DrawLatex(0.95,0.955,'#scale[0.7]{35.6 fb^{-1} (13 TeV)}')

        nuisGrs={}
        dy=1.0/float(len(results)+1)
        for v in varList:
            for ir in range(len(results)):
                title,varVals=results[ir]
                if not title in nuisGrs:               
                    nuisGrs[title]=ROOT.TGraphAsymmErrors()
                    nuisGrs[title].SetTitle(title)
                    nuisGrs[title].SetMarkerStyle(20+ir)
                    nuisGrs[title].SetMarkerColor(COLORS[ir])
                    nuisGrs[title].SetLineColor(COLORS[ir])
                    nuisGrs[title].SetLineWidth(2)
                    nuisGrs[title].SetFillStyle(0)

                if not v in varVals: continue
                npts=nuisGrs[title].GetN()
                val,uncLo,uncHi = varVals[v]
                nuisGrs[title].SetPoint(npts,val,npts+dy*(ir+1))
                nuisGrs[title].SetPointError(npts,abs(uncLo),abs(uncHi),0.,0.)
                
        leg=ROOT.TLegend(0.8,0.9,0.95,0.9-0.06*len(args))
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for title in sorted(nuisGrs):
            nuisGrs[title].Draw('p')
            leg.AddEntry(nuisGrs[title],title,'p')
        leg.Draw()

        c.RedrawAxis()
        c.Modified()
        c.Update()

        for ext in ['png','pdf']:
            c.SaveAs('%s/nuisances_%d.%s'%(outdir,ig,ext))

    return


"""
steer the script
"""
def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #parse user input
    parser = OptionParser(
        usage="%prog -o nuisances label1=fitresults_1.root label2=fitresults_2.root",
        epilog="Summarizes the post-fit results"
        )
    parser.add_option("-o",    type="string",       dest="outdir",  default='nuisances',  help="name of the output directory")
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s'%opt.outdir)
    doNuisanceReport(args,opt.outdir)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())