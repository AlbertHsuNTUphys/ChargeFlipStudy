#!/bin/bash

#parse input arguments
ERA=2017
while getopts "o:y:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then
    echo "Nano_steerAnalysis.sh -o <TEST/SEL> [ -y 2017 ] ";
    echo "  TEST         - test locally the code on a single file";
    echo "  SEL          - launches selection jobs to the batch, output will contain summary trees and control plots";
    echo "  HADD         - produce hadd files";
    exit 1;
fi

#configuration parameters

site=/eos/user/m/melu/TTC_Nanov8_new/
testtag=MC13TeV_2017_DY50toInf_mlm_v2
#testfile=/cms/store/user/tihsu/ttc_nanoNtuple/MC/2017/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_2017/210915_071855/0000/tree_1.root 
testfile=/eos/user/m/melu/TTC_Nanov8_new/DY.root
#testfile=/eos/user/t/tihsu/TTC/MC/2017/DY.root
githash=ttcnanoNtuple
inDir=/cms/store/user/tihsu/ttc_nanoNtuple
outDir=/eos/user/t/tihsu/output_ntuple_2
json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/ExYukawa/input_2017_nano_more.json



#run the operation required

case $WHAT in

    TEST )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${testfile}  --tag ${testtag} \
            -o MC_testsel_${ERA}.root --genWeights genweights_${githash}.root \
            -q local -m RunChargeFlip_nanoaod --era era${ERA} 
        ;;

    SEL )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runeosAnalysis.py \
            --site   ${site} \
            --inDir  ${inDir} \
            --outDir ${outDir} --json ${json}
        ;;

    HADD )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/generate_hadd_file.py \
            --json ${json}
        ;;

esac            
