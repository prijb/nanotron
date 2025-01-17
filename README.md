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
cd nanotron/ && git pull origin Parking && cd ..
```

Build the code
```
scram build clean
scram build
```

List of commands for producing processed Run 2 data trees (common config but with year=2018 or 2018UL option)
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/00496A25-08B6-FB4E-9681-D5FF4E1BE81F.root year=2018 isData=True 

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2430000/0009F6DE-C68E-1E48-8F6D-5603C0B9F215.root year=2018UL isData=True 

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/user/jleonhol/samples/MINIAODSIM/scenarioB1_mpi_4_mA_1p33_ctau_10/MINIAODSIM/231218_121254/0000/miniAOD_1.root year=2018UL isData=False 
```

List of commands for producing processed Run 3 data trees 
```
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2023C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v4/000/367/770/00000/1d32e33a-afa7-45ba-b1c1-a68b9be0409f.root year=2023 isData=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/user/jleonhol/samples/MINIAODSIM/scenarioA_mpi_4_mA_1p33_ctau_10_2022/MINIAODSIM/231122_112812/0000/miniAOD_1.root year=2022 isData=False
```

Test the custom NanoAOD tree (TBD)
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2023
