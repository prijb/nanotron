from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'BParking2023CNew'
config.General.workArea = 'BParking2023CNew'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2023', 'isData=True']

config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v4/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 50
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/BParking2023CNew'
config.Data.publication = True
config.Data.outputDatasetTag = '2023C'

config.Site.storageSite = 'T2_UK_London_IC'