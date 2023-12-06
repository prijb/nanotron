# Main configuration File for custom nanotron "enhanced" nanoAOD ntuples
# 
# CMSSW references:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile
# http://nuclear.ucdavis.edu/~jrobles/CMS/others/SummerStudent2009-1.pdf
#
# m.mieskolainen@imperial.ac.uk, 2023


import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from PhysicsTools.NanoAOD.common_cff import *

from Configuration.StandardSequences.Eras import eras
from RecoBTag.Configuration.RecoBTag_cff import *
from Configuration.AlCa.GlobalTag import GlobalTag


# ------------------------------------------------------------------------
# Processing options

options = VarParsing ('analysis')

options.register(
    'isData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is data"
)

options.register(
    'addSignalLHE',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "adds LHE weights of signal samples"
)

options.register(
    'year',
    '2018',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "add year file"
)

options.register(
    'test',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "running test"
)

options.register(
    'output',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "output filename"
)

options.parseArguments() 

if options.year == '2016':
    process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
elif options.year == '2017':
    process = cms.Process('NANO',eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv2)
elif options.year == '2018' or options.year == '2018D':
    process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)
elif options.year == '2022':
    process = cms.Process('NANO',eras.Run3,eras.run3_nanoAOD_122)
elif options.year == '2023':
    process = cms.Process('NANO',eras.Run3,eras.run3_nanoAOD_124)
else:
    process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)

print("Selected year: ", options.year)

# ------------------------------------------------------------------------
# Import of standard configurations

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

if options.isData:
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    dataTier = cms.untracked.string('NANOAOD')
else:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
    process.load('Configuration.StandardSequences.MagneticField_cff')
    dataTier = cms.untracked.string('NANOAODSIM')

# ------------------------------------------------------------------------
# More options
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

files = {
    'test': {
        "mc": "/store/user/kjpena/miniAODv3_08Feb2020/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_MINIAODSIM/200209_212810/0000/output_1.root",
        #"mc": "https://github.com/LLPDNNX/test-files/raw/master/miniaod/Moriond17_aug2018_miniAODv3_HNL.root",
        },
    '2016': {
        #"mc":"root://xrootd.grid.hep.ph.ic.ac.uk//store/mc/RunIISummer16MiniAODv3/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/40000/0A4EAAB1-9223-E911-B512-A4BF01283A8B.root",
        #"mc":"root://xrootd.grid.hep.ph.ic.ac.uk//store/mc/RunIISummer16MiniAODv3/GluGluToHHTo2B2Tau_node_SM_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/50000/AEF4E98A-C672-E911-9FD1-AC1F6BAC7D18.root",
        #"mc":"root://xrootd.grid.hep.ph.ic.ac.uk//store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/00000/025CA8F7-7C08-E911-8165-0242AC1C0501.root",
        "mc": [
            "root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/miniaod16v3_200929/LLPGun/miniaod16v3_200929/201007_212837/0000/GUN2016_%i.root"%(i) for i in range(1,10)
        ],
        #"mc": "root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/miniaod16v3_200517/HNL_dirac_all_ctau1p0e00_massHNL6p0_Vall6p496e-03/miniaod16v3_200517/200517_004822/0000/HNL2016_140.root", #"root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_tau_Dirac_massiveAndCKM_LO/heavyNeutrino_1.root",
        #"mc": "root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-10_V-0.00112249721603_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_76.root",
        "data": "/store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/00000/16924A85-4D8C-E811-A51C-A4BF01013F29.root",
    },
    '2017': {
        "mc": "root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/testgun17v2/LLPGun/testgun17v2/200724_143441/0000/GUN2017_1.root",
        #"mc": "root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-8_V-0.00214242852856_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root",
        "data": "/store/data/Run2017E/SingleMuon/MINIAOD/31Mar2018-v1/00000/A6325FCE-1C39-E811-BB22-0CC47A745298.root"
    },
    '2018': {
        "mc":"root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/miniaod18_200625/HNL_dirac_all_ctau1p0e-01_massHNL10p0_Vall5p262e-03/miniaod18_200625/200709_103117/0000/HNL2018_283.root",
        #"mc": "root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Autumn18/displaced/HeavyNeutrino_lljj_M-8_V-0.00214242852856_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root",
        "data": "/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FF47BB90-FC1A-CC44-A635-2B8B8C64AA39.root"
    },
    '2018D': {
        "data": "/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FF47BB90-FC1A-CC44-A635-2B8B8C64AA39.root"
    }
}

if len(options.inputFiles) > 0:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles)
    )
else:
    process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring(files[options.year]['data'] if options.isData else files[options.year]['mc'])
        fileNames = cms.untracked.vstring()
    )

# ------------------------------------------------------------------------
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test102X nevts:10000'),
    name       = cms.untracked.string('Applications'),
    version    = cms.untracked.string('$Revision: 1.19 $')
)

# ------------------------------------------------------------------------
# Output definition

output_filename = "nano.root" if not options.output else options.output

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = dataTier,
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
               ['llpnanoAOD_step'] if options.isData else ['llpnanoAOD_step'] # ['llpnanoAOD_step_mu','llpnanoAOD_step_ele'] ~ boolean OR (union) between 'mu' and 'ele' paths
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string(output_filename),
    #outputCommands = process.NANOAODSIMEventContent.outputCommands+cms.untracked.vstring(
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',

        'drop *_caloMetTable_*_*',
        'drop *_saJetTable_*_*',
        'drop *_saTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_fatJetMCTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_subjetMCTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genSubJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',

        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',

        #'drop *_rivetMetTable_*_*',
    )
)

# ------------------------------------------------------------------------
## Output file

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

if options.year == "test":
    options.year = "2016"

if options.isData:
    if options.year == '2016':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v13', '')
    elif options.year == '2017':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v13', '')
    elif options.year == '2018':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v13', '')
    elif options.year == '2018D':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v16', '')
    elif options.year == '2022':
        process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10', '')
    elif options.year == '2023':
        process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_v1', '')
    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None')
else:
    if options.year == '2016':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mcRun2_asymptotic_v8', '')
    elif options.year == '2017':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v8', '')
    elif options.year == '2018' or options.year == '2018D':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21', '')
    elif options.year == '2022':
        process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2022_realistic_postEE_v1', '')
    elif options.year == '2023':
        process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_postBPix_v1', '')
    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

# ------------------------------------------------------------------------
# Custom collections

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#Loading Scouting to PAT sequences
process.load('nanotron.Scouting.scoutingelectron_cff')
process.load('nanotron.Scouting.scoutingjet_cff')
process.load('nanotron.Scouting.scoutingmuon_cff')
process.load('nanotron.Scouting.scoutingvertices_cff')
#process.load('nanotron.Scouting.scoutingpf_cff')


process.muonVerticesTable = cms.EDProducer("MuonVertexProducer",
    srcMuon = cms.InputTag("run3ScoutingMuonToPatMuon"),
    pvSrc   = cms.InputTag("run3ScoutingVertices", "pvs"),
    svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    ptMin   = cms.double(0.8),
    svName  = cms.string("muonSV"),
)

process.fourmuonVerticesTable = cms.EDProducer("FourMuonVertexProducer",
    srcMuon = cms.InputTag("run3ScoutingMuonToPatMuon"),
    pvSrc   = cms.InputTag("run3ScoutingVertices", "pvs"),
    svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    ptMin   = cms.double(0.8),
    svName  = cms.string("fourmuonSV"),
)


process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *  ", # this is the default
    "keep++ pdgId = {Z0}",
    "drop pdgId = {Z0} & status = 2"
    )
)

process.finalGenParticles = cms.EDProducer("GenParticlePruner",
    select = cms.vstring(
        'drop *',
        'keep++ abs(pdgId) == 15 & (pt > 15 ||  isPromptDecayed() )',
        'keep+ abs(pdgId) == 15 ',
        '+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())',
        '+keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15',
        'drop abs(pdgId)= 2212 && abs(pz) > 1000',
        'keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)',
        'keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16',
        'keep status == 3 || (status > 20 && status < 30)',
        'keep isHardProcess() ||  fromHardProcessDecayed()  || fromHardProcessFinalState() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())',
        'keep  (status > 70 && status < 80 && pt > 15) ',
        'keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 37 ',
        'keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)'
    ),
    src = cms.InputTag("prunedGenParticles")
)

# ------------------------------------------------------------------------

# ========================================================================
# ** DATA SEQUENCE **
# ========================================================================

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits


if options.isData:

    # Main
    process.llpnanoAOD_step = cms.Path(
        process.gtStage2Digis + process.l1bits +
        process.electronSequence + process.muonSequence + process.jetSequence
        + process.vertexSequence + process.muonVerticesTable + process.fourmuonVerticesTable
    )

    """
    process.llpnanoAOD_step_ele = cms.Path(
        process.electronFilterSequence+
        process.nanoSequence+
        process.adaptedVertexing+

        process.pfOnionTagInfos+
        process.nanoTable
    )
    """

    # B-parking additions
    # # # # # # # # # # # # # # # # # process.llpnanoAOD_step_mu += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence
    #process.llpnanoAOD_step_ele += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence
    #process.llpnanoAOD_step_mu += process.metadata

# ========================================================================
# ** MC SEQUENCE **
# ========================================================================

else:
    # Main
    process.llpnanoAOD_step = cms.Path(
        process.gtStage2Digis + process.l1bits +
        process.electronSequence + process.muonSequence + process.jetSequence
        + process.vertexSequence + process.muonVerticesTable + process.fourmuonVerticesTable
        #+ process.vertexSequence + process.muonVerticesTable + process.particleSequence
    )
    # # # process.llpnanoAOD_step = cms.Path(
        # # # process.nanoSequenceMC+
        # # # process.adaptedVertexing+
        # # # process.pfOnionTagInfos+
        # # # process.displacedGenVertexSequence+

        # # # process.MCGenDecayInfo+
        # # # process.MCLabels+

        # # # process.nanoTable+
        # # # process.nanoGenTable
    # # # )

    # B-parking additions
    # # # # process.llpnanoAOD_step += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence
    #process.llpnanoAOD_step += process.muonBParkMC # Not used currently

    # LHE
    # # # if options.addSignalLHE:
        # # # process.llpnanoAOD_step += process.lheWeightsTable

    process.finalGenParticlesTask = cms.Task(process.prunedGenParticles, process.finalGenParticles)
    process.mc_path = cms.Path(process.finalGenParticlesTask, process.genParticleTablesTask)

process.endjob_step           = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# ------------------------------------------------------------------------
# Sequence


if options.isData:
#    process.schedule = cms.Schedule(process.llpnanoAOD_step_mu, process.llpnanoAOD_step_ele, process.endjob_step, process.NANOAODSIMoutput_step)
    process.schedule = cms.Schedule(process.llpnanoAOD_step, process.endjob_step, process.NANOAODSIMoutput_step)
else:
    process.schedule = cms.Schedule(process.llpnanoAOD_step, process.mc_path, process.endjob_step, process.NANOAODSIMoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# ------------------------------------------------------------------------
# Remove unneeded modules

modulesToRemove = [
    'jetCorrFactorsAK8',
    'updatedJetsAK8',
    'finalJetsAK8',
    'tightJetIdAK8',
    'looseJetIdAK8',
    'tightJetIdLepVetoAK8',
    'updatedJetsAK8WithUserData',
    'lepInJetVars',
    'chsForSATkJets',
    'softActivityJets',
    'softActivityJets2',
    'softActivityJets5',
    'softActivityJets10',
    'finalJetsAK8',
    'fatJetTable',
    'fatJetMCTable',
    'subJetTable',
    'subjetMCTable',
    'genSubJetAK8Table',
    'saJetTable',
    'saTable',
    "genJetAK8Table",
    "genJetAK8FlavourAssociation",
    "genJetAK8FlavourTable",
   
    "HTXSCategoryTable",
    "rivetProducerHTXS",
    "genSubJetAK8Table",
    
    # "l1bits",
]

#override final jets

#process.finalJets.addBTagInfo=cms.bool(True)
#process.finalJets.addDiscriminators = cms.bool(True)
#process.finalJets.addTagInfos=cms.bool(True)

# # # # # # # # # # for moduleName in modulesToRemove:
    # # # # # # # # # # if hasattr(process,moduleName):
        # # # # # # # # # # print("removing module:", moduleName)
        # # # # # # # # # # if options.isData:
            # # # # # # # # # # process.nanoSequence.remove(getattr(process,moduleName))
        # # # # # # # # # # else:
            # # # # # # # # # # process.nanoSequenceMC.remove(getattr(process,moduleName))
    # # # # # # # # # # else:
        # # # # # # # # # # print("module for removal not found:", moduleName)

#override final photons (required by object linker) so that ID evaluation is not needed
#process.finalPhotons.cut = cms.string("pt > 5")
#process.finalPhotons.src = cms.InputTag("slimmedPhotons")

process.genParticleTable.variables.vertex_x = Var("vertex().x()", float, doc="vertex x position")
process.genParticleTable.variables.vertex_y = Var("vertex().y()", float, doc="vertex y position")
process.genParticleTable.variables.vertex_z = Var("vertex().z()", float, doc="vertex z position")

'''
process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root'),
    outputCommands = process.NANOAODSIMoutput.outputCommands,
    dropMetaData = cms.untracked.string('ALL'),
)
'''
#process.endpath= cms.EndPath(process.OUT)

# ------------------------------------------------------------------------
# Golden lumisection JSON

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
from os import getenv
# if options.isData:
    # import FWCore.PythonUtilities.LumiList as LumiList
    # goldenjson = "nanotron/json/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    # lumilist = LumiList.LumiList(filename=goldenjson).getCMSSWString().split(',')
    # print("Found json list of lumis to process with {} lumi sections from {}".format(len(lumilist), goldenjson))
    # process.source.lumisToProcess = cms.untracked(cms.VLuminosityBlockRange()+lumilist)

# ------------------------------------------------------------------------
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# print(process.dumpPython())
# End adding early deletion
