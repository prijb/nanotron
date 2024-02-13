from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3Scouting2022F'
config.General.workArea = 'Run3Scouting2022F'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanotron/NANOProducer/test/produceNANO.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['year=2022', 'isData=True', 'mode=Scouting']

config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
#Introduce blocks
#2022F blocks
config.Data.inputBlocks = [
    '/ScoutingPFRun3/Run2022F-v1/RAW#0061d92e-19b2-41e6-abae-fdf47ac81aa9',
    '/ScoutingPFRun3/Run2022F-v1/RAW#00be8cd0-0389-4ee8-9e0c-7e0c39bdc94c',
    '/ScoutingPFRun3/Run2022F-v1/RAW#00c579e5-5ea8-4bb5-a49a-b3984e1f82cc',
    '/ScoutingPFRun3/Run2022F-v1/RAW#0104669c-72ff-4193-9af7-6d0a85965df6',
    '/ScoutingPFRun3/Run2022F-v1/RAW#016b6d1e-63d6-469a-98d8-f2ff91999597',
    '/ScoutingPFRun3/Run2022F-v1/RAW#01a8abf0-7949-4a17-bcd0-01fbf622e7d2',
    '/ScoutingPFRun3/Run2022F-v1/RAW#01d84846-1d38-413c-9938-58c37c4643e0',
    '/ScoutingPFRun3/Run2022F-v1/RAW#02389535-169b-4a4f-b315-22fb91723270',
    '/ScoutingPFRun3/Run2022F-v1/RAW#024668b9-5a22-43c1-a301-846e998a93a4',
    '/ScoutingPFRun3/Run2022F-v1/RAW#02a397cb-252d-4427-8b60-9eed599ed4dd',
    '/ScoutingPFRun3/Run2022F-v1/RAW#03073545-163d-4001-906e-b7397ad99265',
    '/ScoutingPFRun3/Run2022F-v1/RAW#03e48820-5db8-4311-a91a-e63936f4a33e',
    '/ScoutingPFRun3/Run2022F-v1/RAW#03f6156a-a94b-4440-9e3c-86f3297f0e82',
    '/ScoutingPFRun3/Run2022F-v1/RAW#04a67ecc-6ee9-4a61-a215-e881a0736b2d',
    '/ScoutingPFRun3/Run2022F-v1/RAW#050663fa-172a-4ff1-9df7-74d7920c7073',
    '/ScoutingPFRun3/Run2022F-v1/RAW#055a898b-218b-4e2b-b4b2-5d1718722e68',
    '/ScoutingPFRun3/Run2022F-v1/RAW#05668658-35f6-4c05-8e9c-6643adfdba11',
    '/ScoutingPFRun3/Run2022F-v1/RAW#05897635-682c-4739-90f1-6acdf9bdc04a',
    '/ScoutingPFRun3/Run2022F-v1/RAW#059a5114-e7fe-4670-84e9-e72f1f281494',
    '/ScoutingPFRun3/Run2022F-v1/RAW#06619b0b-ee93-48e2-9589-d77f44e29d77',
    '/ScoutingPFRun3/Run2022F-v1/RAW#067908ea-af42-4548-b1c7-7ea4f9429187',
    '/ScoutingPFRun3/Run2022F-v1/RAW#069cff88-5daa-4dce-a1db-f25485be2adf',
    '/ScoutingPFRun3/Run2022F-v1/RAW#069e0fb4-30ec-4cd6-a27d-b81f1c6322c3',
    '/ScoutingPFRun3/Run2022F-v1/RAW#06ac0336-7325-4462-8f8d-e4a278b8edc5',
    '/ScoutingPFRun3/Run2022F-v1/RAW#06c9712e-530d-46e4-bb10-ebf26d4f26ac',
    '/ScoutingPFRun3/Run2022F-v1/RAW#06e37e4a-b9d0-4257-b1be-b96f925bb165',
    '/ScoutingPFRun3/Run2022F-v1/RAW#06e7f27f-6592-496c-927f-763c95b3a298'    
]
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 1000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/Run3Scouting2022F'
config.Data.publication = True
config.Data.outputDatasetTag = '2022F'
config.Site.storageSite = 'T2_UK_London_IC'