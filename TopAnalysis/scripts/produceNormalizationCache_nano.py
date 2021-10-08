import os
import sys
import optparse
import subprocess
import json
import ROOT

def get_abbre(name,sample_type,year):
    if sample_type == 'MC':
        return name.split('/')[1] + '_' + year
    elif sample_type == 'data':
        return name.split('/')[1] + '_' + name.split('/')[2].split('-')[0]
  
def fill_histogram(name,sample_type,year,era,subdir,inDir,site, fOut):
    if "Data" in subdir: return
    abbre_name = get_abbre(name,sample_type,year)
    first_name = name.split('/')[1]
    CurrentDir = inDir+'/%s/%s/%s/%s/'%(sample_type,year,first_name,abbre_name)
    print('processing %s'%subdir)

    files = None
    run = 1
    while(run):
        r = subprocess.Popen('xrdfs %s ls %s'%(site,CurrentDir), shell=True, stdout=subprocess.PIPE)
        r_return = str(r.stdout.read().strip()).strip('b')
        if(len(r_return)==2): return
        if 'root' in str(r_return):
            run = 0
            files = r_return.splitlines()
        else:
            CurrentDir = str(r_return)

    h = ROOT.TH1F(subdir,subdir,1,0,1)
    h.SetDirectory(0)

    for f in files:
        fIn = ROOT.TFile.Open(site+f)
        if fIn.IsZombie(): continue
        h_genweight = fIn.Get("nEventsGenWeighted")
        nEventGenweight = h_genweight.GetBinContent(1)
        h.Fill(0.5,nEventGenweight)
        fIn.Close()
    fOut.cd()
    h.Write()



if __name__=='__main__':

  #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s','--site', dest='site', help='remote T2 site', default=None, type='string')
    parser.add_option('-i','--inDir', dest='inDir', help='input directory with files', default=None, type='string')
    parser.add_option('-o','--outFile', dest='outFile', help='output genweights file', default=None, type='string')
    parser.add_option('-f','--json', dest='json', help='json file', default=None, type='string')
    (args,opt) = parser.parse_args()

    print("start producing %s from %s/%s..."%(args.outFile,args.site,args.inDir))
 
    with open(args.json, "r") as f:
        jsons = json.load(f)
        f.close()

    fOut = ROOT.TFile(args.outFile, "RECREATE")
    fOut.cd()
  
    for dataset in jsons:
        fill_histogram(dataset['name'], dataset['type'], str(dataset['year']), dataset['era'], dataset['dir'], args.inDir,args.site,fOut)

    fOut.Close()


