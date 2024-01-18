from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3ScoutingDQCD'
config.General.workArea = 'Run3ScoutingDQCD'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO_scouting.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2022', 'isData=False']

config.Data.inputDataset = '/scenarioA_mpi_4_mA_1p33_ctau_10_2022/jleonhol-MINIAODSIM-dc87d0d68558bef75241a5a1b5633313/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 92
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/Run3ScoutingDQCD'
config.Data.publication = True
config.Data.outputDatasetTag = 'scenarioA_mpi_4_mA_1p33_ctau_10_2022'

config.Site.storageSite = 'T2_UK_London_IC'