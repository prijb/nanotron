from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'BsToJPsiPhi'
config.General.workArea = 'BsToJPsiPhi'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO_AOD.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2022', 'isData=False']

config.Data.inputDataset = '/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 200
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/BsToJPsiPhi'
config.Data.publication = True
config.Data.outputDatasetTag = '2022'

config.Site.storageSite = 'T2_UK_London_IC'