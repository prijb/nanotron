from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'scenarioA_mpi_4_mA_1p33_ctau_10'
config.General.workArea = 'Run3BParkingSignals'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2022preEE', 'isData=False']

config.Data.inputDataset = '/scenarioA_mpi_4_mA_1p33_ctau_10/tafoyava-MINIAODSIM_2022_SecondTry-79e94dc7acbc79867a99aebaabda58b5/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 10000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/samples/Parking/Run3/Nanotronv14'
config.Data.publication = True
config.Data.outputDatasetTag = '2022preEE'

config.Site.storageSite = 'T2_UK_London_IC'