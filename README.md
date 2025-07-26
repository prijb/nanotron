# nanotron

![build tests](https://travis-ci.org/mieskolainen/nanotron.svg?branch=master)

Version 0.02

Based on a fusion of code from:

https://github.com/LLPDNNX/LLPReco

https://github.com/DiElectronX/BParkingNANO

Processes MINIAOD for Run 2 and 3

</br>

Initialize the environment
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Initialize GRID
```
source /cvmfs/grid.cern.ch/alma9-ui-current/etc/profile.d/setup-alma9-test.sh
voms-proxy-init --rfc -voms cms --valid 192:00
```

Initialize CMSSW
```
scram list CMSSW

cmsrel CMSSW_14_0_18
cd CMSSW_14_0_18/src

cmsenv
```

Pull the repo (Version which also has configs)
```
git clone git@github.com:prijb/nanotron.git nanotron 
```

Update the repo
```
cd nanotron/ && git pull origin Parking && cd ..
```

Build the code
```
scram build clean
scram build
```

This version of Nanotron is designed for producing NanoAODv14 ntuples and is incompatible with 2022/2023 data processed using CMSSW 12_4_X (aka Prompt or reprocessed in Dec2022) due to Tau ID changes. The compatible formats are the Sep23 campaigns for 2022, 2023 and the PromptReco for 2024

List of example commands for producing processed Run 3 data trees 
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023-v1/60000/00807cd1-61e0-4a30-917a-f2957e4365b0.root year=2022postEE isData=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2023C/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/0064eaee-25a0-4cf2-a1c2-e33d2d0c43d1.root year=2023preBPix isData=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/000a11f5-055a-4402-b35b-a1e28cabaed7.root year=2023postBPix isData=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2024F/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/382/209/00000/7d2d6f4f-c0b3-4eb5-b351-d1cbb159f308.root year=2024 isData=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2024F/ParkingSingleMuon0/MINIAOD/PromptReco-v1/000/382/213/00000/1787fa93-2b5d-4e2e-bf22-00242c0317a4.root year=2024 isData=True
```

List of example commands for producing processed Run 3 MC trees
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/mc/Run3Summer22MiniAODv4/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/validDigi_130X_mcRun3_2022_realistic_v5-v4/2540000/3b3036b2-c738-4b32-b2f2-b3ba86ad0026.root year=2022preEE isData=False

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/user/jleonhol/samples/MINIAODSIM/scenarioA_mpi_4_mA_1p33_ctau_10_2022/MINIAODSIM/231122_112812/0000/miniAOD_1.root year=2022postEE isData=False
```

Test the custom NanoAOD tree (TBD)
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2023
j.leon-holgado@imperial.ac.uk, 2024
prijith.babu-pradeep18@imperial.ac.uk, 2024