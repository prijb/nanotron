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

Initialize CMSSW

```
scram list CMSSW

cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src

cmsenv
```

Pull the repo
```
git clone git@github.com:mieskolainen/nanotron.git nanotron
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
cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod18.root year=2018 test=True

cmsRun nanotron/NANOProducer/test/produceNANO.py inputFiles=/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/00496A25-08B6-FB4E-9681-D5FF4E1BE81F.root year=2018 test=True isData=True
```

Test the custom NanoAOD tree
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2023
