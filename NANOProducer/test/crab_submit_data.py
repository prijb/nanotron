from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '2022F_updated'
config.General.workArea = 'Run3BParkingData'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2022postEE', 'isData=True']

config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2022F-22Sep2023-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 1
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/Data/Parking/Run3/Nanotronv14'
config.Data.publication = True
config.Data.outputDatasetTag = '2022F_updated'

config.Site.storageSite = 'T2_UK_London_IC'