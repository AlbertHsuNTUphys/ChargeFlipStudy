# ChargeFlip Study

## Download
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_16
cd CMSSW_10_6_16/src
cmsenv
git cms-init
cd $CMSSW_BASE/src
```

## Control Plot
The main code is src/ChargeFlip_nanoaod.cc. To run it, you need to first change the output directory in test/analysis/ExYukawa/Nano_steerAnalysis.sh.
Then to test the code,
```
cmsenv
scram b -j 8
sh test/analysis/ExYukawa/Nano_steerAnalysis.sh -o TEST
```
To run it on condor,
```
cmsenv
scram b -j 8
sh test/analysis/ExYukawa/Nano_steerAnalysis.sh -o SEL
```
To generate the control plot
```
sh plot_and_copy.sh <input_dir> <output_dir>
```


