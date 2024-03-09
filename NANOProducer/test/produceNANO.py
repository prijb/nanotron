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
    'mode',
    'Offline',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Offline or Scouting"
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
    elif options.year == '2018' or options.year == '2018D' or options.year == "2018UL":
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
    elif options.year == '2018' or options.year == '2018D' or options.year == "2018UL":
        #process = cms.Process('NANO',eras.Run2_2018)
        process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_106Xv2)
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
    elif options.year == '2018UL':
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
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
    elif options.year == '2018' or options.year == '2018D' or options.year == "2018UL":
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21', '')
    elif options.year == '2022':
        process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2022_realistic_postEE_v1', '')
    elif options.year == '2023':
        process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_postBPix_v1', '')
    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

# ------------------------------------------------------------------------
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#Loading Scouting to PAT sequences
process.load('nanotron.Scouting.scoutingelectron_cff')
process.load('nanotron.Scouting.scoutingjet_cff')
process.load('nanotron.Scouting.scoutingmuon_cff')
process.load('nanotron.Scouting.scoutingvertices_cff')

#General configs (base + nanotron + muonSV)
#Configs from nano_cff and nanotron loaded only for "offline" 
if options.mode == 'Offline':
    updateJetCollection(
        process,
        labelName      = "OnionTag",
        jetSource      = cms.InputTag('updatedJets'),
        jetCorrections = jetCorrectionsAK4PFchs,
        pfCandidates   = cms.InputTag('packedPFCandidates'),
        pvSource       = cms.InputTag("offlineSlimmedPrimaryVertices"),
        #svSource = cms.InputTag('adaptedSlimmedSecondaryVertices'), 
        svSource       = cms.InputTag('slimmedSecondaryVertices'),
        muSource       = cms.InputTag('slimmedMuons'),
        elSource       = cms.InputTag('slimmedElectrons'),
        btagInfos = [
            # 'pfImpactParameterTagInfos',
            # 'pfInclusiveSecondaryVertexFinderTagInfos',
            # 'pfDeepCSVTagInfos',
            'pfDeepFlavourTagInfos'
        ],
        btagDiscriminators = ['pfDeepFlavourJetTags:probb', 'pfDeepFlavourJetTags:probbb'],
        explicitJTA = False,
    )

    process.pfOnionTagInfos = cms.EDProducer("OnionInfoProducer",
        jets                       = cms.InputTag("updatedJets"),
        muonSrc                    = cms.InputTag("slimmedMuons"),
        electronSrc                = cms.InputTag("slimmedElectrons"),
        shallow_tag_infos          = cms.InputTag('pfDeepCSVTagInfosOnionTag'),
        # shallow_tag_infos          = cms.InputTag('pfDeepFlavourTagInfosOnionTag'), # not usable, dows not produce a shallow_tag_infos
        vertices                   = cms.InputTag('offlineSlimmedPrimaryVertices'),
        secondary_vertices_adapted = cms.InputTag("adaptedSlimmedSecondaryVertices"),
        secondary_vertices         = cms.InputTag("slimmedSecondaryVertices")
    )

    process.nanoTable = cms.EDProducer("NANOProducer",
        srcJets   = cms.InputTag("updatedJets"),
        srcTags   = cms.InputTag("pfOnionTagInfos")
    )

    process.nanoGenTable = cms.EDProducer("NANOGenProducer",
        srcJets   = cms.InputTag("updatedJets"),
        srcLabels = cms.InputTag("MCLabels"),
        srcTags   = cms.InputTag("pfOnionTagInfos")
    )

    process.load('nanotron.NANOProducer.GenDisplacedVertices_cff')

    process.MCGenDecayInfo = cms.EDProducer(
        "MCGenDecayInfoProducer",
        src = cms.InputTag("genParticlesMerged"),
        decays = cms.PSet(
            #dark QCD
            dark_vector = cms.PSet(
                llpId = cms.int32(4900113),
                daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
            ),
            dark_photon = cms.PSet(
                llpId = cms.int32(999999),
                daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
            ),
            #hnl -> qql
            hnl_dirac = cms.PSet(
                llpId = cms.int32(9990012),
                daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
            ),
            hnl_majorana = cms.PSet(
                llpId = cms.int32(9900012),
                daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
            ),
            #gluino -> qq chi0
            split = cms.PSet(
                llpId = cms.int32(1000021),
                daughterIds = cms.vint32([1,2,3,4,5])
            ),
            #gluino -> g gravitino
            gmsb = cms.PSet(
                llpId = cms.int32(1000021),
                daughterIds = cms.vint32([21])
            ),
            #stop -> bl
            rpv = cms.PSet(
                llpId = cms.int32(1000006),
                daughterIds = cms.vint32([5,11,13,15])
            ),
            #H->SS->bbbb
            hss = cms.PSet(
                llpId = cms.int32(9000006),
                daughterIds = cms.vint32([5])
            ),
        )
    )

    process.MCLabels = cms.EDProducer(
        "MCLabelProducer",
        srcVertices  = cms.InputTag("displacedGenVertices"),
        srcJets      = cms.InputTag("finalJets"),
        srcDecayInfo = cms.InputTag("MCGenDecayInfo"),
    )

    process.lheWeightsTable = cms.EDProducer(
        "LHEWeightsProducer",
        lheInfo      = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
        weightGroups = cms.PSet()
    )

    # Particle gun parameters
    process.lheWeightsTable.weightGroups.gun_ctau    = cms.vstring(['ctau'])
    process.lheWeightsTable.weightGroups.gun_llpmass = cms.vstring(['llpmass'])

    # ------------------------------------------------------------------------
    # Coupling reweighting

    process.lheWeightsTable.weightGroups.coupling    = cms.vstring()
    for i in range(1,68):
        process.lheWeightsTable.weightGroups.coupling.append("rwgt_%i"%(i))


    # ------------------------------------------------------------------------
    # Parton Density Functions

    # PDF NNPDF3.1 NNLO hessian
    process.lheWeightsTable.weightGroups.nnpdfhessian = cms.vstring()
    for i in range(1048,1151):
        process.lheWeightsTable.weightGroups.nnpdfhessian.append("%i"%(i))

    # PDF NNPDF3.1 NNLO replicas
    process.lheWeightsTable.weightGroups.nnpdfreplica = cms.vstring()
    for i in range(1151,1252):
        process.lheWeightsTable.weightGroups.nnpdfreplica.append("%i"%(i))

    # Scale weights
    for scaleSet in [
        ['murNominal_mufNominal',range(1001,1006)],
        ['murUp_mufNominal',     range(1006,1011)],
        ['murDown_mufNominal',   range(1011,1016)],
        ['murNominal_mufUp',     range(1016,1021)],
        ['murUp_mufUp',          range(1021,1026)],
        ['murDown_mufUp',        range(1026,1031)],
        ['murNominal_mufDown',   range(1031,1036)],
        ['murUp_mufDown',        range(1036,1041)],
        ['murDown_mufDown',      range(1041,1046)],
        ['emission',             range(1046,1048)],
        ]:
        setattr(process.lheWeightsTable.weightGroups,scaleSet[0],cms.vstring())
        for i in scaleSet[1]:
            getattr(process.lheWeightsTable.weightGroups,scaleSet[0]).append("%i"%(i))


    process.load('RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff')
    process.load('nanotron.NANOProducer.adaptedSV_cff')

    # ------------------------------------------------------------------------
    # Customisation function from PhysicsTools.NanoAOD.nano_cff
    # Note: customizeCommon is essential since it breaks the nanotronSequence process afterwards without it!

    # from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData, nanoAOD_customizeMC
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon

    if options.isData:
        # process = nanoAOD_customizeData(process)
        process = nanoAOD_customizeCommon(process)
    else:
        # process = nanoAOD_customizeMC(process)
        process = nanoAOD_customizeCommon(process)

    #All the nanotron additions which may or may not work
    process.nanotronSequence = cms.Sequence(process.adaptedVertexing
                                + process.pfOnionTagInfos
                                + process.nanoTable
    ) 
    
    process.nanotronMCSequence = cms.Sequence(process.displacedGenVertexSequence
                                + process.MCGenDecayInfo
                                + process.MCLabels
                                + process.nanoGenTable
    )
    
    process.muonVerticesTable = cms.EDProducer("MuonVertexProducer",
        srcMuon = cms.InputTag("linkedObjects", "muons"),
        pvSrc   = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
        dlenMin = cms.double(0),
        dlenSigMin = cms.double(0),
        ptMin   = cms.double(0.8),
        svName  = cms.string("muonSV"),
    )

    process.fourmuonVerticesTable = cms.EDProducer("FourMuonVertexProducer",
        srcMuon = cms.InputTag("linkedObjects", "muons"),
        #srcMuon = cms.InputTag("finalMuons"),
        pvSrc   = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
        dlenMin = cms.double(0),
        dlenSigMin = cms.double(0),
        ptMin   = cms.double(0.8),
        svName  = cms.string("fourmuonSV"),
    )

#Adding the gen table and muon matching tasks manually from genparticles_cff and muons_cff
if options.mode == 'Scouting':
    process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
    process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

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
    
    process.scoutingSequence = cms.Sequence(process.gtStage2Digis + process.l1bits + process.electronSequence 
        + process.muonSequence + process.jetSequence + process.vertexSequence)

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

#Trigger matching (only for offline)
if options.mode == 'Offline':
    Path=["HLT_Mu7_IP4","HLT_Mu8_IP6","HLT_Mu8_IP5","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","HLT_Mu9_IP6","HLT_Mu9_IP5","HLT_Mu9_IP4","HLT_Mu10p5_IP3p5","HLT_Mu12_IP6"]

    if options.year in ['2022', '2023']:
        Path = [
            'HLT_Dimuon0_Jpsi3p5_Muon2',
            'HLT_Dimuon0_Jpsi_L1_4R_0er1p5R',
            'HLT_Dimuon0_Jpsi_L1_NoOS',
            'HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R',
            'HLT_Dimuon0_Jpsi_NoVertexing_NoOS',
            'HLT_Dimuon0_Jpsi_NoVertexing',
            'HLT_Dimuon0_Jpsi',
            'HLT_Dimuon0_LowMass_L1_0er1p5R',
            'HLT_Dimuon0_LowMass_L1_0er1p5',
            'HLT_Dimuon0_LowMass_L1_4R',
            'HLT_Dimuon0_LowMass_L1_4',
            'HLT_Dimuon0_LowMass_L1_TM530',
            'HLT_Dimuon0_LowMass',
            'HLT_Dimuon0_Upsilon_L1_4p5',
            'HLT_Dimuon0_Upsilon_L1_4p5er2p0M',
            'HLT_Dimuon0_Upsilon_L1_4p5er2p0',
            'HLT_Dimuon0_Upsilon_Muon_NoL1Mass',
            'HLT_Dimuon0_Upsilon_NoVertexing',
            'HLT_Dimuon10_Upsilon_y1p4',
            'HLT_Dimuon12_Upsilon_y1p4',
            'HLT_Dimuon14_Phi_Barrel_Seagulls',
            'HLT_Dimuon14_PsiPrime_noCorrL1',
            'HLT_Dimuon14_PsiPrime',
            'HLT_Dimuon18_PsiPrime_noCorrL1',
            'HLT_Dimuon18_PsiPrime',
            'HLT_Dimuon24_Phi_noCorrL1',
            'HLT_Dimuon24_Upsilon_noCorrL1',
            'HLT_Dimuon25_Jpsi_noCorrL1',
            'HLT_Dimuon25_Jpsi',
            'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05',
            'HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon',
            'HLT_DoubleMu3_TkMu_DsTau3Mu',
            'HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass',
            'HLT_DoubleMu3_Trk_Tau3mu',
            'HLT_DoubleMu4_3_Bs',
            'HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG',
            'HLT_DoubleMu4_3_Jpsi',
            'HLT_DoubleMu4_3_LowMass',
            'HLT_DoubleMu4_3_Photon4_BsToMMG',
            'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
            'HLT_DoubleMu4_JpsiTrk_Bc',
            'HLT_DoubleMu4_Jpsi_Displaced',
            'HLT_DoubleMu4_Jpsi_NoVertexing',
            'HLT_DoubleMu4_LowMass_Displaced',
            'HLT_DoubleMu4_MuMuTrk_Displaced',
            'HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL',
            'HLT_Mu25_TkMu0_Phi',
            'HLT_Mu30_TkMu0_Psi',
            'HLT_Mu30_TkMu0_Upsilon',
            'HLT_Mu4_L1DoubleMu',
            'HLT_Mu7p5_L2Mu2_Jpsi',
            'HLT_Mu7p5_L2Mu2_Upsilon',
            'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1',
            'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15',
            'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1',
            'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15',
            'HLT_Trimuon5_3p5_2_Upsilon_Muon',
            'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon'
        ]
    
    process.muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                #muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD    
                                muonCollection = cms.InputTag("linkedObjects", "muons"), #same collection as in muonSV                                                     
                                bits = cms.InputTag("TriggerResults","","HLT"),
                                prescales = cms.InputTag("patTrigger"),
                                objects = cms.InputTag("slimmedPatTrigger"),
                                vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),

                                ##for the output trigger matched collection
                                maxdR_matching = cms.double(0.1),

                                ## for the output selected collection (tag + all compatible in dZ)
                                #filterMuon is redundant since this version doesn't skip muons
                                filterMuon = cms.bool(True), 
                                dzForCleaning_wrtTrgMuon = cms.double(1.0),

                                ptMin = cms.double(0.5),
                                absEtaMax = cms.double(2.4),
                                # keeps only muons with at soft Quality flag
                                softMuonsOnly = cms.bool(False),
                                HLTPaths=cms.vstring(Path)#, ### comma to the softMuonsOnly
    #				 L1seeds=cms.vstring(Seed)
                                )
    
    process.muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("muonTrgSelector:SelectedMuons"),
        cut = cms.string(""), #we should not filter on cross linked collections
        name = cms.string("MuonBPark"),
        doc  = cms.string("slimmedMuons for BPark after basic selection"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(CandVars,
            ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
            dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
            dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
            dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
            dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
            vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
            vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
            vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
            ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
            sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
    #        segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
    #        nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
            #nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
    #        pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
            pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
            pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
    #        tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
            isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
            isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
            isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
            mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
    #        mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
            tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
            softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
    #        softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
    #        highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
            pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
            tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
    #        mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
    #        miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
    #        multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
            triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
    #        inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
            isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco muon is also triggering"),#########################################################3,
    #        toWhichHLTisMatched = Var("userInt('ToWhichHLTisMatched')",int,doc="To which HLT muons is the reco muon matched, -1 for probe" ),
            matched_dr = Var("userFloat('DR')",float,doc="dr with the matched triggering muon" ),
            matched_dpt = Var("userFloat('DPT')",float,doc="dpt/pt with the matched triggering muon" ),        #comma
            passCuts = Var("userInt('passCuts')",bool,doc="Does the muon pass the pt and eta cuts, and is matched in dZ with the trigger muon?"),
            looseId = Var("userInt('looseId')",int,doc="reco muon is Loose"),
            # fired_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5')",int,doc="reco muon fired this trigger"),
            # fired_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6')",int,doc="reco muon fired this trigger")#,
        ),
    )

    for p in Path:
        setattr(process.muonBParkTable.variables, "fired_%s" % p, Var("userInt('%s')" % p, int, doc="reco muon fired this trigger"))

    process.muonTriggerMatchedTable = process.muonBParkTable.clone(
        src = cms.InputTag("muonTrgSelector:trgMuons"),
        name = cms.string("TriggerMuon"),
        doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
        variables = cms.PSet(CandVars,
            vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
            vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
            vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
            ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
            dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
            dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
            dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
            dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
            ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
            sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10)
    #        trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
        )
    )

    trigobjpaths = ['HLT_DoubleMu4_3_LowMass', 'HLT_DoubleMu4_LowMass_Displaced', 'HLT_Dimuon10_Upsilon_y1p4']

    #Introduce trigger muons without any L1 selections
    process.triggerMuonTable = cms.EDProducer("TriggerObjectProducer",
        bits = cms.InputTag("TriggerResults","","HLT"),
        objects = cms.InputTag("slimmedPatTrigger"),
        name= cms.string("TriggerObject"),
        ptMin = cms.double(0.5),
        objId = cms.int32(83),
        HLTPaths = cms.vstring(trigobjpaths)
        
    )

    from  PhysicsTools.NanoAOD.triggerObjects_cff import *

    process.triggerObjectBParkTable = cms.EDProducer("TriggerObjectTableBParkProducer",
        name= cms.string("TrigObjBPark"),
        src = cms.InputTag("unpackedPatTrigger"),
        l1Muon = cms.InputTag("gmtStage2Digis","Muon"),
        selections = cms.VPSet(
            cms.PSet(
                name = cms.string("Muon"),
                id = cms.int32(13),
                sel = cms.string("type(83) && pt > 5 && coll('hltIterL3MuonCandidates')"), 
                l1seed = cms.string("type(-81)"), l1deltaR = cms.double(0.5),
                l2seed = cms.string("type(83) && coll('hltL2MuonCandidates')"),  l2deltaR = cms.double(0.3),
                qualityBits = cms.string("filter('hltL3fL1s*Park*')"), qualityBitsDoc = cms.string("1 = Muon filters for BPH parking"),
            ),
        )
    )
    
    process.muonBParkSequence  = cms.Sequence(process.muonTrgSelector)# * process.countTrgMuons)
    process.muonBParkTables    = cms.Sequence(process.muonBParkTable)
    process.muonTriggerMatchedTables = cms.Sequence(process.muonTriggerMatchedTable)   ####
    process.triggerObjectBParkTables = cms.Sequence(unpackedPatTrigger + process.triggerObjectBParkTable + process.triggerMuonTable)

# ------------------------------------------------------------------------
#Order of sequences is dependent on dataset and whether it's data or MC
if options.mode == 'Offline':
    if options.isData:
        process.llpnanoAOD_step = cms.Path(
            process.nanoSequence
            + process.nanotronSequence
        )
    else:
        process.llpnanoAOD_step = cms.Path(
            process.nanoSequenceMC
            + process.nanotronSequence
            #+ process.nanotronMCSequence    #Debugging
        )

        # LHE
        if options.addSignalLHE:
            process.llpnanoAOD_step += process.lheWeightsTable
    
    process.llpnanoAOD_step += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence

if options.mode == 'Scouting':
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

# Remove unneeded modules in offline

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
    "genSubJetAK8Table"
]

#Remove more modules (for 2018?)
#Check how much we can get in common for data vs MC
if (options.year == '2018' or options.year == "2018UL"):
    #Delete processes common to data and MC
    del process.lhcInfoTable

    #Modify the collections
    #Simplifying tau cut because some of the MVAs don't work
    run2_nanoAOD_106Xv2.toModify(
        process.finalTaus,
        cut = cms.string("pt > 18")
    )
    #Gets rid of low pT electrons in cross linker
    run2_nanoAOD_106Xv2.toModify(
        process.linkedObjects, lowPtElectrons=None
    )

    if options.year == '2018':
        #Using slimmed electrons instead of final electrons in cross linker since in case adding MVA user data breaks things -> Make sure to remove the electron table if you do this
        run2_nanoAOD_106Xv2.toModify(
            process.linkedObjects, electrons=cms.InputTag("slimmedElectrons")
        )
        del process.electronTable
    
    if options.year == '2018UL':
        #Modifying electron MVA assignment
        run2_nanoAOD_106Xv2.toModify(
            process.electronMVATTH.variables,
            LepGood_pt = cms.string("pt"),
            LepGood_eta = cms.string("eta"),
            LepGood_jetNDauChargedMVASel = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0"),
            LepGood_miniRelIsoCharged = None,
            LepGood_miniRelIsoNeutral = None,
            LepGood_jetPtRelv2 = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0"),
            LepGood_jetDF = cms.string("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0"),
            LepGood_jetPtRatio = None,
            LepGood_dxy = cms.string("log(abs(dB('PV2D')))"),
            LepGood_sip3d = cms.string("abs(dB('PV3D')/edB('PV3D'))"),
            LepGood_dz = cms.string("log(abs(dB('PVDZ')))"),
            LepGood_mvaFall17V2noIso = None
        )

    #Removing problematic jet variables
    del process.jetTable.variables.hfadjacentEtaStripsSize
    del process.jetTable.variables.hfcentralEtaStripSize
    del process.jetTable.variables.hfsigmaEtaEta
    del process.jetTable.variables.hfsigmaPhiPhi
    #Remove proton tables
    del process.singleRPTable
    del process.multiRPTable
    del process.protonTable
    #Remove low pT electrons
    del process.lowPtElectronTable
    #Could probably preserve the electron table if I removed the MVA variables one by one or fixed the MVA assignment
    del process.isoTrackTable

    if options.isData:
        print("Removing data specific modules")

    else:
        print("Removing MC specific modules")    
        del process.genWeightsTable
        del process.genTable
        del process.particleLevel
        del process.lheWeightsTable
        del process.electronMCTable
        del process.lowPtElectronMCTable
        del process.genParticles2HepMC
        del process.genParticles2HepMCHiggsVtx
        del process.tautagger
        del process.rivetLeptonTable
        del process.rivetPhotonTable
        del process.rivetMetTable
        del process.HTXSCategoryTable
        del process.rivetProducerHTXS

if options.mode == 'Offline':
    for moduleName in modulesToRemove:
        if hasattr(process,moduleName):
            print("removing module:", moduleName)
            if options.isData:
                process.nanoSequence.remove(getattr(process,moduleName))
            else:
                process.nanoSequenceMC.remove(getattr(process,moduleName))
        else:
            print("module for removal not found:", moduleName)

process.genParticleTable.variables.vertex_x = Var("vertex().x()", float, doc="vertex x position")
process.genParticleTable.variables.vertex_y = Var("vertex().y()", float, doc="vertex y position")
process.genParticleTable.variables.vertex_z = Var("vertex().z()", float, doc="vertex z position")

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