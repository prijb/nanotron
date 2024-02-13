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

List of commands for producing processed Run 3 data trees (common config with option for Offline and Scouting datasets)
Use format=AOD for files that are in AODSIM
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2023C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v4/000/367/770/00000/1d32e33a-afa7-45ba-b1c1-a68b9be0409f.root year=2023 isData=True mode=Offline

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/user/jleonhol/samples/MINIAODSIM/scenarioA_mpi_4_mA_1p33_ctau_10_2022/MINIAODSIM/231122_112812/0000/miniAOD_1.root year=2022 isData=False mode=Offline

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=file:/vols/cms/pb4918/StoreNTuple/Scouting/2022FSample.root year=2022 isData=True mode=Scouting

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/user/jleonhol/samples/MINIAODSIM/scenarioA_mpi_4_mA_1p33_ctau_10_2022/MINIAODSIM/231122_112812/0000/miniAOD_1.root year=2022 isData=False mode=Scouting

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/mc/Run3Summer22EEDRPremix/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/0b63e2c8-9381-445f-b963-7237608d249e.root year=2022 isData=False mode=Scouting format=AOD
```

Test the custom NanoAOD tree (TBD)
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2023
