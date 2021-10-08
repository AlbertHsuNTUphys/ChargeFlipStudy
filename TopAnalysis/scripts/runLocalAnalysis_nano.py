import os
import sys
import optparse
import subprocess
import json

def get_abbre(name,sample_type,year):
    if sample_type == 'MC':
        return name.split('/')[1] + '_' + year
    elif sample_type == 'data':
        return name.split('/')[1] + '_' + name.split('/')[2].split('-')[0]
  
def prepare_shell(name,sample_type,year,era,subdir,inDir,outDir,site, FarmDir,condor,total_shell):
    abbre_name = get_abbre(name,sample_type,year)
    first_name = name.split('/')[1]
    CurrentDir = inDir+'/%s/%s/%s/%s/'%(sample_type,year,first_name,abbre_name)
    print('processing %s'%subdir)
    run = 1
    files = None
    while(run):
        r = subprocess.Popen('xrdfs %s ls %s'%(site,CurrentDir), shell=True, stdout=subprocess.PIPE)
        r_return = str(r.stdout.read().strip()).strip('b')
        if(len(r_return)==2): return
        if 'root' in str(r_return):
            run = 0
            files = r_return.splitlines()
#            print(files)
        else:
            CurrentDir = str(r_return)

    cmsswBase=os.environ['CMSSW_BASE']

    for f in files:
#        print(f)
        f = f.strip('\'')
        name = f.split('/')[-1]
        index = (name.strip('.root')).strip('tree_')
        shell_name = subdir+'_'+index+'.sh'
        fname = subdir+'_'+index+'.root'
        with open('%s/%s'%(FarmDir,shell_name),'w') as shell:
            shell.write('#!/bin/bash\n')
            shell.write('WORKDIR=`pwd`\n')
            shell.write('cd %s\n'%cmsswBase)
            shell.write('eval `scram r -sh`\n')
            shell.write('cd ${WORKDIR}\n')
            shell.write('python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py '%cmsswBase)
            shell.write('-i %s/%s '%(site,f))
            shell.write('-o ${WORKDIR}/%s '%fname)
            shell.write('--era %s/src/TopLJets2015/TopAnalysis/data/era2017 '%cmsswBase)
            shell.write('--genWeights genweights_ttcnanoNtuple.root ')
            shell.write('--tag %s '%subdir)
            shell.write('--method RunChargeFlip_nanoaod\n')
            shell.write('mv -v ${WORKDIR}/%s %s/%s'%(fname,outDir,fname))
        total_shell.write('sh %s\n'%shell_name)
        condor.write('cfgFile=%s\n'%shell_name)
        condor.write('queue 1\n')





if __name__=='__main__':

  #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s','--site', dest='site', help='remote T2 site', default=None, type='string')
    parser.add_option('-i','--inDir', dest='inDir', help='input directory with files', default=None, type='string')
    parser.add_option('-o','--outDir', dest='outDir', help='output directory with files', default=None, type='string')
    parser.add_option('-f','--json', dest='json', help='json file', default=None, type='string')
    (args,opt) = parser.parse_args()
    print("start merging...")
    FarmDir = os.environ['CMSSW_BASE']+"/NanoAnalysis_FARM"
#    subprocess.Popen(f'rm -rf {FarmDir}', shell=True,stdout=subprocess.PIPE)
    os.system('mkdir -p %s'%FarmDir)
    os.system('rm %s/*'%FarmDir)
    os.system('mkdir -p %s'%args.outDir)
    os.system('rm %s/*.root'%args.outDir)
    proxy = subprocess.Popen('voms-proxy-info -path', shell=True, stdout=subprocess.PIPE)
    proxy_path = ((str(proxy.stdout.read()).strip('b')).splitlines())[0]
    os.system('cp %s %s/x509up'%(proxy_path,FarmDir))

    with open(args.json, "r") as f:
        jsons = json.load(f)
        f.close()

    condor=open('%s/condor.sub'%FarmDir,'w')
    condor.write('use_x509userproxy = true\n')
    condor.write('x509userproxy = %s/x509up\n'%FarmDir)
    condor.write('output = %s/job_common.out\n'%FarmDir)
    condor.write('error  = %s/job_common.err\n'%FarmDir)
    condor.write('log    = %s/job_common.log\n'%FarmDir)
    condor.write('executable = %s/$(cfgFile)\n'%FarmDir)
    condor.write('requirements = (OpSysAndVer =?= "CentOS7")\n')
    condor.write('+JobFlavour = "tomorrow"\n')

    total_shell=open('%s/total_shell.sh'%FarmDir,'w')
  
    for dataset in jsons:
        prepare_shell(dataset['name'], dataset['type'], str(dataset['year']), dataset['era'], dataset['dir'], args.inDir, args.outDir,args.site,FarmDir,condor,total_shell)

    condor.close()
    total_shell.close()
    os.system('condor_submit %s/condor.sub'%FarmDir)


