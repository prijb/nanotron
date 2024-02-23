from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ZeroBias'
config.General.workArea = 'ZeroBias'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2023', 'isData=True', 'mode=Offline']

config.Data.inputDataset = '/ZeroBias/Run2023D-22Sep2023_v1-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 346
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/ZeroBias'
config.Data.publication = True
config.Data.outputDatasetTag = '2023D'

config.Site.storageSite = 'T2_UK_London_IC'