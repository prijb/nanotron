# nanotron

![build tests](https://travis-ci.org/mieskolainen/nanotron.svg?branch=master)

Version 0.02

Based on a fusion of code from:

https://github.com/LLPDNNX/LLPReco

https://github.com/DiElectronX/BParkingNANO

</br>

Initialize the environment
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Initialize GRID
```
source /vols/grid/cms/setup.sh
voms-proxy-init --voms cms
```

Initialize CMSSW
```
scram list CMSSW

cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src

cmsenv
```

Pull the repo (Version which also has configs)
```
git clone git@github.com:prijb/nanotron.git nanotron 
```

Update the repo
```
cd nanotron/ && git pull origin main && cd ..
```

Build the code
```
scram build clean
scram build
```

Produce a custom NanoAOD tree
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod18.root year=2018 test=True isData=False

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/100000/005594DA-4AA0-3E48-A8C2-46DECDE2E925.root year=2018 test=True isData=False addSignalLHE=False

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/00496A25-08B6-FB4E-9681-D5FF4E1BE81F.root year=2018 test=True isData=True
```

List of commands for producing processed Run 3 data trees 
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2023C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v4/000/367/770/00000/1d32e33a-afa7-45ba-b1c1-a68b9be0409f.root year=2023 isData=True

cmsRun nanotron/NANOProducer/test/produceNANO_scouting.py inputFiles=file:/vols/cms/pb4918/StoreNTuple/Scouting/2022FSample.root year=2022 isData=True

cmsRun nanotron/NANOProducer/test/produceNANO_AOD.py inputFiles=/store/mc/Run3Summer22EEDRPremix/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/0b63e2c8-9381-445f-b963-7237608d249e.root year=2022 isData=False
```

Test the custom NanoAOD tree (TBD)
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2023
