# nanotron

![build tests](https://travis-ci.org/mieskolainen/nanotron.svg?branch=master)

Version 0.01

Hard fork based on https://github.com/LLPDNNX/LLPReco


Initialize the environment
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Pull the repo
```
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src
cmsenv
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
```

Test the custom NanoAOD tree
```
python nanotron/test/check_tree.py
```


m.mieskolainen@imperial.ac.uk, 2022
