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
    False,
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
    'format',
    'MINI',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MINI or AOD"
)

options.parseArguments() 

#Caution: Choose eras based on MC/Data
#Example: eras.Run3 works for MC, but eras.Run3,eras.run3_nanoAOD_122 for Data

if options.isData:
    if options.year == '2016':
        process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
    elif options.year == '2017':
        process = cms.Process('NANO',eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv2)
    elif options.year == '2018' or options.year == '2018D':
        process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_106Xv2)
    elif options.year == '2022':
        process = cms.Process('NANO',eras.Run3,eras.run3_nanoAOD_122)
    elif (options.year == '2023'):
        process = cms.Process('NANO',eras.Run3,eras.run3_nanoAOD_124)
    else:
        process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
else:
    if options.year == '2016':
        process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
    elif options.year == '2017':
        process = cms.Process('NANO',eras.Run2_2017)
    elif options.year == '2018' or options.year == '2018D':
        process = cms.Process('NANO',eras.Run2_2018)
    elif options.year == '2022':
        process = cms.Process('NANO',eras.Run3)
    elif (options.year == '2023'):
        process = cms.Process('NANO',eras.Run3)
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
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#The else here makes it compatible with CRAB submissions
if len(options.inputFiles) > 0:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
        bypassVersionCheck = cms.untracked.bool(True)
    )
else:
    process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring(files[options.year]['data'] if options.isData else files[options.year]['mc'])
        fileNames = cms.untracked.vstring(),
        bypassVersionCheck = cms.untracked.bool(True)
    )

# ------------------------------------------------------------------------
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test102X nevts:10000'),
    name       = cms.untracked.string('Applications'),
    version    = cms.untracked.string('$Revision: 1.19 $')
)

# ------------------------------------------------------------------------
# Output definition (NanoAODOutputModule) 

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
               ['llpnanoAOD_step'] #Can make this data/MC dependent
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string('nano.root'),
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
## Output file (What does this even do? Isn't this pointless?)

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

# ------------------------------------------------------------------------
## GlobalTag definition

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
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.l1trig_cff import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#Loading Scouting to PAT sequences
process.load('nanotron.Scouting.scoutingelectron_cff')
process.load('nanotron.Scouting.scoutingjet_cff')
process.load('nanotron.Scouting.scoutingmuon_cff')
process.load('nanotron.Scouting.scoutingvertices_cff')
process.load('nanotron.Scouting.scoutingglobal_cff')

#Get GT digis
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

#Extract L1 muons from the GT (Note: NOT GMT)
process.l1MuScoutingTable = l1MuTable.clone(src = cms.InputTag("gtStage2Digis","Muon"))
process.l1MuScoutingTable.variables = cms.PSet(l1MuonReducedVars)
process.l1MuonSequence = cms.Sequence(process.l1MuScoutingTable)

#Different source from normal configs (probably should just clone for simplicity)
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

#Genmatching sequence 
process.muonsMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = (process.ScoutingMuonTable).src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticles"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)

process.muonMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = (process.ScoutingMuonTable).src,
    mcMap   = cms.InputTag("muonsMCMatchForTable"),
    objName = (process.ScoutingMuonTable).name,
    objType = (process.ScoutingMuonTable).name, #cms.string("Muon"),
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

process.scoutingSequence = cms.Sequence(process.l1MuonSequence + process.electronSequence 
    + process.muonSequence + process.jetSequence + process.vertexSequence + process.globalSequence)

#MC sequence depends on whether the input format is MINIAODSIM or AODSIM
if options.format == 'AOD':
    #PAT conversion and Pruning necessary. Note: separate sequence since tasks and sequences can't be mixed in the same line 
    process.load("Configuration.StandardSequences.PAT_cff")
    process.patSequence = cms.Sequence(process.patTask)
    process.mcSequence = cms.Sequence(process.patSequence + process.prunedGenParticles + process.finalGenParticles + process.genParticleTable + process.muonsMCMatchForTable + process.muonMCTable)

else:    
    process.mcSequence = cms.Sequence(process.finalGenParticles + process.genParticleTable + process.muonsMCMatchForTable + process.muonMCTable)

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

process.muonVertexSequence = cms.Sequence(process.muonVerticesTable + process.fourmuonVerticesTable)

# ------------------------------------------------------------------------
#Order of sequences is dependent on dataset and whether it's data or MC
if options.isData:
    process.llpnanoAOD_step = cms.Path(process.gtStage2Digis 
        + process.l1bits 
        + process.scoutingSequence
        + process.muonVertexSequence
    )
else:
    process.llpnanoAOD_step = cms.Path(process.gtStage2Digis 
        + process.l1bits 
        + process.scoutingSequence
        + process.muonVertexSequence
        + process.mcSequence
    )

    # LHE
    if options.addSignalLHE:
        process.llpnanoAOD_step += process.lheWeightsTable


process.endjob_step           = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# ------------------------------------------------------------------------
#Final sequence
process.schedule = cms.Schedule(process.llpnanoAOD_step, process.endjob_step, process.NANOAODSIMoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.genParticleTable.variables.vertex_x = Var("vertex().x()", float, doc="vertex x position")
process.genParticleTable.variables.vertex_y = Var("vertex().y()", float, doc="vertex y position")
process.genParticleTable.variables.vertex_z = Var("vertex().z()", float, doc="vertex z position")

# ------------------------------------------------------------------------
# Golden lumisection JSON

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
from os import getenv

# ------------------------------------------------------------------------
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# print(process.dumpPython())