function run_test()
{
    cd ~
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc700 || return 1
    scramv1 project CMSSW CMSSW_10_2_22 || return 1
    cd CMSSW_10_2_22/src || return 1
    eval `scram runtime -sh` || return 1
    rsync -r /scripts/* nanotron  || return 1
    scram b || return 1
    pip install --user -r nanotron/requirements.txt || return 1
    
    cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod18.root year=2018 test=True || return 1
    python nanotron/test/check_nan.py || return 1
}

run_test
