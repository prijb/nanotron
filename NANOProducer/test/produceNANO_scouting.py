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
        fileNames = cms.untracked.vstring(files[options.year]['data'] if options.isData else files[options.year]['mc'])
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


process.run3ScoutingJetToPatJet = cms.EDProducer("Run3ScoutingJetToPatJetProducer",
    jetSource=cms.InputTag("hltScoutingPFPacker"))

process.run3ScoutingEleToPatEle = cms.EDProducer("Run3ScoutingEleToPatEleProducer",
    eleSource=cms.InputTag("hltScoutingEgammaPacker"))

process.run3ScoutingPhotonToPatPhoton = cms.EDProducer("Run3ScoutingPhotonToPatPhotonProducer",
    photonSource=cms.InputTag("hltScoutingEgammaPacker"))

process.run3ScoutingMuonRecoTrack = cms.EDProducer("Run3ScoutingMuonRecoTrackProducer",
    muonSource=cms.InputTag("hltScoutingMuonPacker"))

process.run3ScoutingMuonToPatMuon = cms.EDProducer("Run3ScoutingMuonToPatMuonProducer",
    muonSource=cms.InputTag("hltScoutingMuonPacker"),
    trackSource=cms.InputTag("run3ScoutingMuonRecoTrack")
)

process.run3ScoutingVertices = cms.EDProducer("Run3ScoutingVtxToVtxProducer",
    pvSource=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    svSource=cms.InputTag("hltScoutingMuonPacker", "displacedVtx"),
    muonSource=cms.InputTag("hltScoutingMuonPacker")
)

# process.linkedObjects = cms.EDProducer("PATObjectCrossLinker",
    # boostedTaus = cms.InputTag(""),
    # electrons = cms.InputTag("run3ScoutingEleToPatEle"),
    # jets = cms.InputTag("run3ScoutingJetToPatJet", "jets"),
    # lowPtElectrons = cms.InputTag(""),
    # muons = cms.InputTag("run3ScoutingMuonToPatMuon"),
    # photons = cms.InputTag("run3ScoutingPhotonToPatPhoton"),
    # taus = cms.InputTag("run3ScoutingJetToPatJet", "taus"),
    # vertices = cms.InputTag("run3ScoutingVertices", "svs")
# )


process.jetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string('slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15)'),
    extension = cms.bool(False),
    # externalVariables = cms.PSet(
        # bRegCorr = cms.PSet(
            # doc = cms.string('pt correction for b-jet energy regression'),
            # precision = cms.int32(10),
            # src = cms.InputTag("bjetNN","corr"),
            # type = cms.string('float')
        # ),
        # bRegRes = cms.PSet(
            # doc = cms.string('res on pt corrected with b-jet regression'),
            # precision = cms.int32(6),
            # src = cms.InputTag("bjetNN","res"),
            # type = cms.string('float')
        # ),
        # cRegCorr = cms.PSet(
            # doc = cms.string('pt correction for c-jet energy regression'),
            # precision = cms.int32(10),
            # src = cms.InputTag("cjetNN","corr"),
            # type = cms.string('float')
        # ),
        # cRegRes = cms.PSet(
            # doc = cms.string('res on pt corrected with c-jet regression'),
            # precision = cms.int32(6),
            # src = cms.InputTag("cjetNN","res"),
            # type = cms.string('float')
        # )
    # ),
    maxLen = cms.optional.uint32,
    mightGet = cms.optional.untracked.vstring,
    name = cms.string('Jet'),
    singleton = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    # src = cms.InputTag("linkedObjects", "jets"),
    src = cms.InputTag("run3ScoutingJetToPatJet", "jets"),
    variables = cms.PSet(
        area = cms.PSet(
            doc = cms.string('jet catchment area, for JECs'),
            expr = cms.string('1'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagDeepFlavB = cms.PSet(
            doc = cms.string('DeepJet b+bb+lepb tag discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagDeepFlavCvB = cms.PSet(
            doc = cms.string('DeepJet c vs b+bb+lepb discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagDeepFlavCvL = cms.PSet(
            doc = cms.string('DeepJet c vs uds+g discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagDeepFlavQG = cms.PSet(
            doc = cms.string('DeepJet g vs uds discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagRobustParTAK4B = cms.PSet(
            doc = cms.string('RobustParTAK4 b+bb+lepb tag discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagRobustParTAK4CvB = cms.PSet(
            doc = cms.string('RobustParTAK4 c vs b+bb+lepb discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagRobustParTAK4CvL = cms.PSet(
            doc = cms.string('RobustParTAK4 c vs uds+g discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        btagRobustParTAK4QG = cms.PSet(
            doc = cms.string('RobustParTAK4 g vs uds discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        chEmEF = cms.PSet(
            doc = cms.string('charged Electromagnetic Energy Fraction'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        chFPV0EF = cms.PSet(
            doc = cms.string('charged fromPV==0 Energy Fraction (energy excluded from CHS jets). Previously called betastar.'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        chHEF = cms.PSet(
            doc = cms.string('charged Hadron Energy Fraction'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        electronIdx1 = cms.PSet(
            doc = cms.string('index of first matching electron'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        electronIdx2 = cms.PSet(
            doc = cms.string('index of second matching electron'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        eta = cms.PSet(
            doc = cms.string('eta'),
            expr = cms.string('eta'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        hfadjacentEtaStripsSize = cms.PSet(
            doc = cms.string('eta size of the strips next to the central tower strip in HF (noise discriminating variable) '),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        hfcentralEtaStripSize = cms.PSet(
            doc = cms.string('eta size of the central tower strip in HF (noise discriminating variable) '),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        hfsigmaEtaEta = cms.PSet(
            doc = cms.string('sigmaEtaEta for HF jets (noise discriminating variable)'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        hfsigmaPhiPhi = cms.PSet(
            doc = cms.string('sigmaPhiPhi for HF jets (noise discriminating variable)'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        jetId = cms.PSet(
            doc = cms.string('Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        mass = cms.PSet(
            doc = cms.string('mass'),
            expr = cms.string('mass'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        muEF = cms.PSet(
            doc = cms.string('muon Energy Fraction'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        muonIdx1 = cms.PSet(
            doc = cms.string('index of first matching muon'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        muonIdx2 = cms.PSet(
            doc = cms.string('index of second matching muon'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        muonSubtrFactor = cms.PSet(
            doc = cms.string('1-(muon-subtracted raw pt)/(raw pt)'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        nConstituents = cms.PSet(
            doc = cms.string('Number of particles in the jet'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        nElectrons = cms.PSet(
            doc = cms.string('number of electrons in the jet'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        nMuons = cms.PSet(
            doc = cms.string('number of muons in the jet'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        nSVs = cms.PSet(
            doc = cms.string('number of secondary vertices in the jet'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        neEmEF = cms.PSet(
            doc = cms.string('neutral Electromagnetic Energy Fraction'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        neHEF = cms.PSet(
            doc = cms.string('neutral Hadron Energy Fraction'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        phi = cms.PSet(
            doc = cms.string('phi'),
            expr = cms.string('phi'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        pt = cms.PSet(
            doc = cms.string('pt'),
            expr = cms.string('pt'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        puId = cms.PSet(
            doc = cms.string('Pileup ID flags with 106X (2018) training'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        puIdDisc = cms.PSet(
            doc = cms.string('Pileup ID discriminant with 106X (2018) training'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        qgl = cms.PSet(
            doc = cms.string('Quark vs Gluon likelihood discriminator'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        rawFactor = cms.PSet(
            doc = cms.string('1 - Factor to get back to raw pT'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        svIdx1 = cms.PSet(
            doc = cms.string('index of first matching secondary vertex'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        svIdx2 = cms.PSet(
            doc = cms.string('index of second matching secondary vertex'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        )
    )
)

process.electronTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string('slimmedElectrons after basic selection (pt > 5 )'),
    extension = cms.bool(False),
    # externalVariables = cms.PSet(
        # fsrPhotonIdx = cms.PSet(
            # doc = cms.string('Index of the lowest-dR/ET2 among associated FSR photons'),
            # precision = cms.int32(-1),
            # src = cms.InputTag("leptonFSRphotons","eleFsrIndex"),
            # type = cms.string('int16')
        # ),
        # mvaTTH = cms.PSet(
            # doc = cms.string('TTH MVA lepton ID score'),
            # precision = cms.int32(14),
            # src = cms.InputTag("electronMVATTH"),
            # type = cms.string('float')
        # )
    # ),
    maxLen = cms.optional.uint32,
    mightGet = cms.optional.untracked.vstring,
    name = cms.string('Electron'),
    singleton = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    # src = cms.InputTag("linkedObjectsEle","electrons"),
    src = cms.InputTag("run3ScoutingEleToPatEle"),
    variables = cms.PSet(
        charge = cms.PSet(
            doc = cms.string('electric charge'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        convVeto = cms.PSet(
            doc = cms.string('pass conversion veto'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        cutBased = cms.PSet(
            doc = cms.string('cut-based ID RunIII Winter22 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        cutBased_Fall17V2 = cms.PSet(
            doc = cms.string('cut-based ID Fall17V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        cutBased_HEEP = cms.PSet(
            doc = cms.string('cut-based HEEP ID'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        deltaEtaSC = cms.PSet(
            doc = cms.string('delta eta (SC,ele) with sign'),
            expr = cms.string('1'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dr03EcalRecHitSumEt = cms.PSet(
            doc = cms.string('Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        dr03HcalDepth1TowerSumEt = cms.PSet(
            doc = cms.string('Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        dr03TkSumPt = cms.PSet(
            doc = cms.string('Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        dr03TkSumPtHEEP = cms.PSet(
            doc = cms.string('Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV used in HEEP ID'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        dxy = cms.PSet(
            doc = cms.string('dxy (with sign) wrt first PV, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dxyErr = cms.PSet(
            doc = cms.string('dxy uncertainty, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        dz = cms.PSet(
            doc = cms.string('dz (with sign) wrt first PV, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dzErr = cms.PSet(
            doc = cms.string('dz uncertainty, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        eInvMinusPInv = cms.PSet(
            doc = cms.string('1/E_SC - 1/p_trk'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        energyErr = cms.PSet(
            doc = cms.string('energy error of the cluster-track combination'),
            expr = cms.string("1"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        eta = cms.PSet(
            doc = cms.string('eta'),
            expr = cms.string('eta'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        hoe = cms.PSet(
            doc = cms.string('H over E'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        ip3d = cms.PSet(
            doc = cms.string('3D impact parameter wrt first PV, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        isPFcand = cms.PSet(
            doc = cms.string('electron is PF candidate'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        jetIdx = cms.PSet(
            doc = cms.string('index of the associated jet (-1 if none)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        jetNDauCharged = cms.PSet(
            doc = cms.string('number of charged daughters of the closest jet'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        jetPtRelv2 = cms.PSet(
            doc = cms.string('Relative momentum of the lepton with respect to the closest jet after subtracting the lepton'),
            expr = cms.string("1"),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        jetRelIso = cms.PSet(
            doc = cms.string('Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)'),
            expr = cms.string("1"),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        jetRelIso_Fall17V2 = cms.PSet(
            doc = cms.string('Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)'),
            expr = cms.string("1"),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        lostHits = cms.PSet(
            doc = cms.string('number of missing inner hits'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        mass = cms.PSet(
            doc = cms.string('mass'),
            expr = cms.string('mass'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        miniPFRelIso_all = cms.PSet(
            doc = cms.string('mini PF relative isolation, total (with scaled rho*EA PU Winter22V1 corrections)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        miniPFRelIso_all_Fall17V2 = cms.PSet(
            doc = cms.string('mini PF relative isolation, total (with scaled rho*EA PU Fall17V2 corrections)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        miniPFRelIso_chg = cms.PSet(
            doc = cms.string('mini PF relative isolation, charged component'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        miniPFRelIso_chg_Fall17V2 = cms.PSet(
            doc = cms.string('mini PF relative isolation, charged component'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaHZZIso = cms.PSet(
            doc = cms.string('HZZ MVA Iso ID score'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaIso = cms.PSet(
            doc = cms.string('MVA Iso ID score, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaIso_Fall17V2 = cms.PSet(
            doc = cms.string('MVA Iso ID score, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaIso_Fall17V2_WP80 = cms.PSet(
            doc = cms.string('MVA Iso ID WP80, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaIso_Fall17V2_WP90 = cms.PSet(
            doc = cms.string('MVA Iso ID WP90, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaIso_Fall17V2_WPL = cms.PSet(
            doc = cms.string('MVA Iso ID loose WP, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaIso_WP80 = cms.PSet(
            doc = cms.string('MVA Iso ID WP80, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaIso_WP90 = cms.PSet(
            doc = cms.string('MVA Iso ID WP90, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaNoIso = cms.PSet(
            doc = cms.string('MVA noIso ID score, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaNoIso_Fall17V2 = cms.PSet(
            doc = cms.string('MVA noIso ID score, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        mvaNoIso_Fall17V2_WP80 = cms.PSet(
            doc = cms.string('MVA noIso ID WP80, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaNoIso_Fall17V2_WP90 = cms.PSet(
            doc = cms.string('MVA noIso ID WP90, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaNoIso_Fall17V2_WPL = cms.PSet(
            doc = cms.string('MVA noIso ID loose WP, Fall17V2'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaNoIso_WP80 = cms.PSet(
            doc = cms.string('MVA noIso ID WP80, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mvaNoIso_WP90 = cms.PSet(
            doc = cms.string('MVA noIso ID WP90, Winter22V1'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        pdgId = cms.PSet(
            doc = cms.string('PDG code assigned by the event reconstruction (not by MC truth)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        pfRelIso03_all = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3, total (with rho*EA PU Winter22V1 corrections)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pfRelIso03_all_Fall17V2 = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3 with 94 EffArea, total (with rho*EA PU corrections)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pfRelIso03_chg = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3, charged component'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pfRelIso03_chg_Fall17V2 = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3 with 94 EffArea, charged component'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        phi = cms.PSet(
            doc = cms.string('phi'),
            expr = cms.string('phi'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        photonIdx = cms.PSet(
            doc = cms.string('index of the first associated photon (-1 if none)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        pt = cms.PSet(
            doc = cms.string('pt'),
            expr = cms.string('pt'),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        r9 = cms.PSet(
            doc = cms.string('R9 of the supercluster, calculated with full 5x5 region'),
            expr = cms.string('1'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        scEtOverPt = cms.PSet(
            doc = cms.string('(supercluster transverse energy)/pt-1'),
            expr = cms.string('1'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        seedGain = cms.PSet(
            doc = cms.string('Gain of the seed crystal'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        seediEtaOriX = cms.PSet(
            doc = cms.string('iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100.'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('int8')
        ),
        seediPhiOriY = cms.PSet(
            doc = cms.string('iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100.'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        sieie = cms.PSet(
            doc = cms.string('sigma_IetaIeta of the supercluster, calculated with full 5x5 region'),
            expr = cms.string('1'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        sip3d = cms.PSet(
            doc = cms.string('3D impact parameter significance wrt first PV, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        svIdx = cms.PSet(
            doc = cms.string('index of matching secondary vertex'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        tightCharge = cms.PSet(
            doc = cms.string('Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        vidNestedWPBitmap = cms.PSet(
            doc = cms.string('VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleEBEECut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEBEECut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        vidNestedWPBitmapHEEP = cms.PSet(
            doc = cms.string('VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleFull5x5SigmaIEtaIEtaWithSatCut,GsfEleFull5x5E2x5OverE5x5WithSatCut,GsfEleHadronicOverEMLinearCut,GsfEleTrkPtIsoCut,GsfEleEmHadD1IsoRhoCut,GsfEleDxyCut,GsfEleMissingHitsCut,GsfEleEcalDrivenCut), 1 bits per cut'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        vidNestedWPBitmap_Fall17V2 = cms.PSet(
            doc = cms.string('VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleEBEECut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEBEECut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int')
        )
    )
)

process.muonTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string("slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed(\'CutBasedIdLoose\') || passed(\'SoftCutBasedId\') || passed(\'SoftMvaId\') || passed(\'CutBasedIdGlobalHighPt\') || passed(\'CutBasedIdTrkHighPt\'))))"),
    extension = cms.bool(False),
    # externalVariables = cms.PSet(
        # fsrPhotonIdx = cms.PSet(
            # doc = cms.string('Index of the lowest-dR/ET2 among associated FSR photons'),
            # precision = cms.int32(-1),
            # src = cms.InputTag("leptonFSRphotons","muFsrIndex"),
            # type = cms.string('int16')
        # ),
        # mvaLowPt = cms.PSet(
            # doc = cms.string('Low pt muon ID score'),
            # precision = cms.int32(14),
            # src = cms.InputTag("muonMVALowPt"),
            # type = cms.string('float')
        # ),
        # mvaTTH = cms.PSet(
            # doc = cms.string('TTH MVA lepton ID score'),
            # precision = cms.int32(14),
            # src = cms.InputTag("muonMVATTH"),
            # type = cms.string('float')
        # )
    # ),
    maxLen = cms.optional.uint32,
    mightGet = cms.optional.untracked.vstring,
    name = cms.string('Muon'),
    singleton = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    src = cms.InputTag("run3ScoutingMuonToPatMuon"),
    variables = cms.PSet(
        charge = cms.PSet(
            doc = cms.string('electric charge'),
            expr = cms.string('charge'),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        dxy = cms.PSet(
            doc = cms.string('dxy (with sign) wrt first PV, in cm'),
            expr = cms.string("dB(\'PV2D\')"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dxyErr = cms.PSet(
            doc = cms.string('dxy uncertainty, in cm'),
            expr = cms.string("edB(\'PV2D\')"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        dxybs = cms.PSet(
            doc = cms.string('dxy (with sign) wrt the beam spot, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dz = cms.PSet(
            doc = cms.string('dz (with sign) wrt first PV, in cm'),
            expr = cms.string("dB(\'PVDZ\')"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        dzErr = cms.PSet(
            doc = cms.string('dz uncertainty, in cm'),
            expr = cms.string("abs(edB(\'PVDZ\'))"),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        eta = cms.PSet(
            doc = cms.string('eta'),
            expr = cms.string('eta'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        highPtId = cms.PSet(
            doc = cms.string('high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        highPurity = cms.PSet(
            doc = cms.string('inner track is high purity'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        inTimeMuon = cms.PSet(
            doc = cms.string('inTimeMuon ID'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        ip3d = cms.PSet(
            doc = cms.string('3D impact parameter wrt first PV, in cm'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        isGlobal = cms.PSet(
            doc = cms.string('muon is global muon'),
            expr = cms.string('isGlobalMuon'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        isPFcand = cms.PSet(
            doc = cms.string('muon is PF candidate'),
            expr = cms.string('isPFMuon'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        isStandalone = cms.PSet(
            doc = cms.string('muon is a standalone muon'),
            expr = cms.string('isStandAloneMuon'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        isTracker = cms.PSet(
            doc = cms.string('muon is tracker muon'),
            expr = cms.string('isTrackerMuon'),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        jetIdx = cms.PSet(
            doc = cms.string('index of the associated jet (-1 if none)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        jetNDauCharged = cms.PSet(
            doc = cms.string('number of charged daughters of the closest jet'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        jetPtRelv2 = cms.PSet(
            doc = cms.string('Relative momentum of the lepton with respect to the closest jet after subtracting the lepton'),
            expr = cms.string("1"),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        jetRelIso = cms.PSet(
            doc = cms.string('Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)'),
            expr = cms.string("1"),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        looseId = cms.PSet(
            doc = cms.string('muon is loose muon'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mass = cms.PSet(
            doc = cms.string('mass'),
            expr = cms.string('mass'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        mediumId = cms.PSet(
            doc = cms.string('cut-based ID, medium WP'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        mediumPromptId = cms.PSet(
            doc = cms.string('cut-based ID, medium prompt WP'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        miniIsoId = cms.PSet(
            doc = cms.string('MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        miniPFRelIso_all = cms.PSet(
            doc = cms.string('mini PF relative isolation, total (with scaled rho*EA PU corrections)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        miniPFRelIso_chg = cms.PSet(
            doc = cms.string('mini PF relative isolation, charged component'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        multiIsoId = cms.PSet(
            doc = cms.string('MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        mvaMuID = cms.PSet(
            doc = cms.string('MVA-based ID score '),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        mvaMuID_WP = cms.PSet(
            doc = cms.string('MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        nStations = cms.PSet(
            doc = cms.string('number of matched stations with default arbitration (segment & track)'),
            # expr = cms.string("userInt('numberofmatchedstations')"),
            # expr = cms.string("numberofmatchedstations"),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        nTrackerLayers = cms.PSet(
            doc = cms.string('number of layers in the tracker'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        pdgId = cms.PSet(
            doc = cms.string('PDG code assigned by the event reconstruction (not by MC truth)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        pfIsoId = cms.PSet(
            doc = cms.string('PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        pfRelIso03_all = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3, total (deltaBeta corrections)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pfRelIso03_chg = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.3, charged component'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pfRelIso04_all = cms.PSet(
            doc = cms.string('PF relative isolation dR=0.4, total (deltaBeta corrections)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        phi = cms.PSet(
            doc = cms.string('phi'),
            expr = cms.string('phi'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        pt = cms.PSet(
            doc = cms.string('pt'),
            expr = cms.string('pt'),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        ptErr = cms.PSet(
            doc = cms.string('ptError of the muon track'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        puppiIsoId = cms.PSet(
            doc = cms.string('PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        segmentComp = cms.PSet(
            doc = cms.string('muon segment compatibility'),
            expr = cms.string('1'),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        sip3d = cms.PSet(
            doc = cms.string('3D impact parameter significance wrt first PV'),
            expr = cms.string("1"),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        softId = cms.PSet(
            doc = cms.string('soft cut-based ID'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        softMva = cms.PSet(
            doc = cms.string('soft MVA ID score'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        softMvaId = cms.PSet(
            doc = cms.string('soft MVA ID'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        svIdx = cms.PSet(
            doc = cms.string('index of matching secondary vertex'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('int16')
        ),
        tightCharge = cms.PSet(
            doc = cms.string('Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)'),
            expr = cms.string('1'),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        tightId = cms.PSet(
            doc = cms.string('cut-based ID, tight WP'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        tkIsoId = cms.PSet(
            doc = cms.string('TkIso ID (1=TkIsoLoose, 2=TkIsoTight)'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        tkRelIso = cms.PSet(
            doc = cms.string('Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        triggerIdLoose = cms.PSet(
            doc = cms.string('TriggerIdLoose ID'),
            expr = cms.string("1"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        tunepRelPt = cms.PSet(
            doc = cms.string('TuneP relative pt, tunePpt/pt'),
            expr = cms.string('1'),
            precision = cms.int32(6),
            type = cms.string('float')
        ),
        vx = cms.PSet(
            doc = cms.string('Vertex x'),
            expr = cms.string('track().vx()'),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        vy = cms.PSet(
            doc = cms.string('Vertex y'),
            expr = cms.string('track().vy()'),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        vz = cms.PSet(
            doc = cms.string('Vertex z'),
            expr = cms.string('track().vz()'),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
    )
)

process.vertexTable = cms.EDProducer("VertexTableProducer",
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(3),
    goodPvCut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
    pvName = cms.string('PV'),
    pvSrc = cms.InputTag("run3ScoutingVertices", "pvs"),
    svCut = cms.string(''),
    svDoc = cms.string('secondary vertices from IVF algorithm'),
    svName = cms.string('SV'),
    svSrc = cms.InputTag("run3ScoutingVertices", "svs")
)

# process.pvTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    # cut = cms.string(''),
    # doc = cms.string(''),
    # extension = cms.bool(True),
    # externalVariables = cms.PSet(

    # ),
    # maxLen = cms.optional.uint32,
    # mightGet = cms.optional.untracked.vstring,
    # name = cms.string('SV'),
    # singleton = cms.bool(False),
    # skipNonExistingSrc = cms.bool(False),
    # src = cms.InputTag("run3ScoutingVertices", "pvs"),
    # variables = cms.PSet(
        # chi2 = cms.PSet(
            # doc = cms.string('reduced chi2, i.e. chi/ndof'),
            # expr = cms.string('vertexNormalizedChi2()'),
            # precision = cms.int32(8),
            # type = cms.string('float')
        # ),
        # eta = cms.PSet(
            # doc = cms.string('eta'),
            # expr = cms.string('eta'),
            # precision = cms.int32(12),
            # type = cms.string('float')
        # ),
        # mass = cms.PSet(
            # doc = cms.string('mass'),
            # expr = cms.string('1'),
            # precision = cms.int32(10),
            # type = cms.string('float')
        # ),
        # ndof = cms.PSet(
            # doc = cms.string('number of degrees of freedom'),
            # expr = cms.string('vertexNdof()'),
            # precision = cms.int32(8),
            # type = cms.string('float')
        # ),
        # ntracks = cms.PSet(
            # doc = cms.string('number of tracks'),
            # expr = cms.string('numberOfDaughters()'),
            # precision = cms.int32(-1),
            # type = cms.string('uint8')
        # ),
        # phi = cms.PSet(
            # doc = cms.string('phi'),
            # expr = cms.string('1'),
            # precision = cms.int32(12),
            # type = cms.string('float')
        # ),
        # pt = cms.PSet(
            # doc = cms.string('pt'),
            # expr = cms.string('1'),
            # precision = cms.int32(10),
            # type = cms.string('float')
        # ),
        # x = cms.PSet(
            # doc = cms.string('primary vertex X position, in cm'),
            # expr = cms.string('position().x()'),
            # precision = cms.int32(10),
            # type = cms.string('float')
        # ),
        # y = cms.PSet(
            # doc = cms.string('primary vertex Y position, in cm'),
            # expr = cms.string('position().y()'),
            # precision = cms.int32(10),
            # type = cms.string('float')
        # ),
        # z = cms.PSet(
            # doc = cms.string('primary vertex Z position, in cm'),
            # expr = cms.string('position().z()'),
            # precision = cms.int32(14),
            # type = cms.string('float')
        # )
    # )
# )

process.svTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string(''),
    extension = cms.bool(True),
    externalVariables = cms.PSet(

    ),
    maxLen = cms.optional.uint32,
    mightGet = cms.optional.untracked.vstring,
    name = cms.string('SV'),
    singleton = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    # src = cms.InputTag("run3ScoutingVertices", "svs"),
    src = cms.InputTag("vertexTable"),
    variables = cms.PSet(
        chi2 = cms.PSet(
            doc = cms.string('reduced chi2, i.e. chi/ndof'),
            expr = cms.string('vertexNormalizedChi2()'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        eta = cms.PSet(
            doc = cms.string('eta'),
            expr = cms.string('eta'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        mass = cms.PSet(
            doc = cms.string('mass'),
            expr = cms.string('mass'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        ndof = cms.PSet(
            doc = cms.string('number of degrees of freedom'),
            expr = cms.string('vertexNdof()'),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        ntracks = cms.PSet(
            doc = cms.string('number of tracks'),
            expr = cms.string('numberOfDaughters()'),
            precision = cms.int32(-1),
            type = cms.string('uint8')
        ),
        phi = cms.PSet(
            doc = cms.string('phi'),
            expr = cms.string('phi'),
            precision = cms.int32(12),
            type = cms.string('float')
        ),
        pt = cms.PSet(
            doc = cms.string('pt'),
            expr = cms.string('pt'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        x = cms.PSet(
            doc = cms.string('primary vertex X position, in cm'),
            expr = cms.string('position().x()'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        y = cms.PSet(
            doc = cms.string('primary vertex Y position, in cm'),
            expr = cms.string('position().y()'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        z = cms.PSet(
            doc = cms.string('primary vertex Z position, in cm'),
            expr = cms.string('position().z()'),
            precision = cms.int32(14),
            type = cms.string('float')
        )
    )
)




# process.TablesTask = cms.Task(process.linkedObjects, process.linkedObjectsEle, process.jetTable, process.electronTable)
process.electronSequence = cms.Sequence((process.run3ScoutingEleToPatEle * process.electronTable) + process.run3ScoutingPhotonToPatPhoton)
process.muonSequence = cms.Sequence(process.run3ScoutingMuonRecoTrack * process.run3ScoutingMuonToPatMuon * process.muonTable)
process.jetSequence = cms.Sequence(process.run3ScoutingJetToPatJet * process.jetTable)
# process.vertexSequence = cms.Sequence((process.run3ScoutingPVtoVertex * process.pvTable) + (process.run3ScoutingSVtoVertex * process.svTable))
process.vertexSequence = cms.Sequence(process.run3ScoutingVertices * process.vertexTable * (process.svTable))


process.muonVerticesTable = cms.EDProducer("MuonVertexProducer",
    srcMuon = cms.InputTag("run3ScoutingMuonToPatMuon"),
    pvSrc   = cms.InputTag("run3ScoutingVertices", "pvs"),
    svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    ptMin   = cms.double(0.8),
    svName  = cms.string("muonSV"),
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

# process.vertexSequence = cms.Sequence((process.run3ScoutingPVtoVertex + process.run3ScoutingSVtoVertex) * process.vertexTable)
# process.TablesTask = process.electronSequence + process.muonSequence + process.jetSequence

# # # # # # # # # # updateJetCollection(
    # # # # # # # # # # process,
    # # # # # # # # # # labelName      = "OnionTag",
    # # # # # # # # # # jetSource      = cms.InputTag('updatedJets'),
    # # # # # # # # # # jetCorrections = jetCorrectionsAK4PFchs,
    # # # # # # # # # # pfCandidates   = cms.InputTag('packedPFCandidates'),
    # # # # # # # # # # pvSource       = cms.InputTag("offlineSlimmedPrimaryVertices"),
    # # # # # # # # # # #svSource = cms.InputTag('adaptedSlimmedSecondaryVertices'), 
    # # # # # # # # # # svSource       = cms.InputTag('slimmedSecondaryVertices'),
    # # # # # # # # # # muSource       = cms.InputTag('slimmedMuons'),
    # # # # # # # # # # elSource       = cms.InputTag('slimmedElectrons'),
    # # # # # # # # # # btagInfos = [
        # # # # # # # # # # # 'pfImpactParameterTagInfos',
        # # # # # # # # # # # 'pfInclusiveSecondaryVertexFinderTagInfos',
        # # # # # # # # # # # 'pfDeepCSVTagInfos',
        # # # # # # # # # # 'pfDeepFlavourTagInfos'
    # # # # # # # # # # ],
    # # # # # # # # # # btagDiscriminators = ['pfDeepFlavourJetTags:probb', 'pfDeepFlavourJetTags:probbb'],
    # # # # # # # # # # explicitJTA = False,
# # # # # # # # # # )


# # # # # # # # # # process.pfOnionTagInfos = cms.EDProducer("OnionInfoProducer",
    # # # # # # # # # # jets                       = cms.InputTag("updatedJets"),
    # # # # # # # # # # muonSrc                    = cms.InputTag("slimmedMuons"),
    # # # # # # # # # # electronSrc                = cms.InputTag("slimmedElectrons"),
    # # # # # # # # # # shallow_tag_infos          = cms.InputTag('pfDeepCSVTagInfosOnionTag'),
    # # # # # # # # # # # shallow_tag_infos          = cms.InputTag('pfDeepFlavourTagInfosOnionTag'), # not usable, dows not produce a shallow_tag_infos
    # # # # # # # # # # vertices                   = cms.InputTag('offlineSlimmedPrimaryVertices'),
    # # # # # # # # # # secondary_vertices_adapted = cms.InputTag("adaptedSlimmedSecondaryVertices"),
    # # # # # # # # # # secondary_vertices         = cms.InputTag("slimmedSecondaryVertices")
# # # # # # # # # # )


# # # # # # # # # # process.nanoTable = cms.EDProducer("NANOProducer",
    # # # # # # # # # # srcJets   = cms.InputTag("updatedJets"),
    # # # # # # # # # # srcTags   = cms.InputTag("pfOnionTagInfos")
# # # # # # # # # # )


# # # # # # # # # # process.nanoGenTable = cms.EDProducer("NANOGenProducer",
    # # # # # # # # # # srcJets   = cms.InputTag("updatedJets"),
    # # # # # # # # # # srcLabels = cms.InputTag("MCLabels"),
    # # # # # # # # # # srcTags   = cms.InputTag("pfOnionTagInfos")
# # # # # # # # # # )


# # # # # # # # # # process.load('nanotron.NANOProducer.GenDisplacedVertices_cff')

# # # # # # # # # # process.MCGenDecayInfo = cms.EDProducer(
    # # # # # # # # # # "MCGenDecayInfoProducer",
    # # # # # # # # # # src = cms.InputTag("genParticlesMerged"),
    # # # # # # # # # # decays = cms.PSet(
        # # # # # # # # # # #dark QCD
        # # # # # # # # # # dark_vector = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(4900113),
            # # # # # # # # # # daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
        # # # # # # # # # # ),
        # # # # # # # # # # dark_photon = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(999999),
            # # # # # # # # # # daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
        # # # # # # # # # # ),
        # # # # # # # # # # #hnl -> qql
        # # # # # # # # # # hnl_dirac = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(9990012),
            # # # # # # # # # # daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
        # # # # # # # # # # ),
        # # # # # # # # # # hnl_majorana = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(9900012),
            # # # # # # # # # # daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
        # # # # # # # # # # ),
        # # # # # # # # # # #gluino -> qq chi0
        # # # # # # # # # # split = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(1000021),
            # # # # # # # # # # daughterIds = cms.vint32([1,2,3,4,5])
        # # # # # # # # # # ),
        # # # # # # # # # # #gluino -> g gravitino
        # # # # # # # # # # gmsb = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(1000021),
            # # # # # # # # # # daughterIds = cms.vint32([21])
        # # # # # # # # # # ),
        # # # # # # # # # # #stop -> bl
        # # # # # # # # # # rpv = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(1000006),
            # # # # # # # # # # daughterIds = cms.vint32([5,11,13,15])
        # # # # # # # # # # ),
        # # # # # # # # # # #H->SS->bbbb
        # # # # # # # # # # hss = cms.PSet(
            # # # # # # # # # # llpId = cms.int32(9000006),
            # # # # # # # # # # daughterIds = cms.vint32([5])
        # # # # # # # # # # ),
    # # # # # # # # # # )
# # # # # # # # # # )

# # # # # # # # # # process.MCLabels = cms.EDProducer(
    # # # # # # # # # # "MCLabelProducer",
    # # # # # # # # # # srcVertices  = cms.InputTag("displacedGenVertices"),
    # # # # # # # # # # srcJets      = cms.InputTag("finalJets"),
    # # # # # # # # # # srcDecayInfo = cms.InputTag("MCGenDecayInfo"),
# # # # # # # # # # )

# # # # # # # # # # process.lheWeightsTable = cms.EDProducer(
    # # # # # # # # # # "LHEWeightsProducer",
    # # # # # # # # # # lheInfo      = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
    # # # # # # # # # # weightGroups = cms.PSet()
# # # # # # # # # # )

# # # # # # # # # # # Particle gun parameters
# # # # # # # # # # process.lheWeightsTable.weightGroups.gun_ctau    = cms.vstring(['ctau'])
# # # # # # # # # # process.lheWeightsTable.weightGroups.gun_llpmass = cms.vstring(['llpmass'])

# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # Coupling reweighting

# # # # # # # # # # process.lheWeightsTable.weightGroups.coupling    = cms.vstring()
# # # # # # # # # # for i in range(1,68):
    # # # # # # # # # # process.lheWeightsTable.weightGroups.coupling.append("rwgt_%i"%(i))


# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # Parton Density Functions

# # # # # # # # # # # PDF NNPDF3.1 NNLO hessian
# # # # # # # # # # process.lheWeightsTable.weightGroups.nnpdfhessian = cms.vstring()
# # # # # # # # # # for i in range(1048,1151):
    # # # # # # # # # # process.lheWeightsTable.weightGroups.nnpdfhessian.append("%i"%(i))

# # # # # # # # # # # PDF NNPDF3.1 NNLO replicas
# # # # # # # # # # process.lheWeightsTable.weightGroups.nnpdfreplica = cms.vstring()
# # # # # # # # # # for i in range(1151,1252):
    # # # # # # # # # # process.lheWeightsTable.weightGroups.nnpdfreplica.append("%i"%(i))

# # # # # # # # # # # Scale weights
# # # # # # # # # # for scaleSet in [
    # # # # # # # # # # ['murNominal_mufNominal',range(1001,1006)],
    # # # # # # # # # # ['murUp_mufNominal',     range(1006,1011)],
    # # # # # # # # # # ['murDown_mufNominal',   range(1011,1016)],
    # # # # # # # # # # ['murNominal_mufUp',     range(1016,1021)],
    # # # # # # # # # # ['murUp_mufUp',          range(1021,1026)],
    # # # # # # # # # # ['murDown_mufUp',        range(1026,1031)],
    # # # # # # # # # # ['murNominal_mufDown',   range(1031,1036)],
    # # # # # # # # # # ['murUp_mufDown',        range(1036,1041)],
    # # # # # # # # # # ['murDown_mufDown',      range(1041,1046)],
    # # # # # # # # # # ['emission',             range(1046,1048)],
    # # # # # # # # # # ]:
    # # # # # # # # # # setattr(process.lheWeightsTable.weightGroups,scaleSet[0],cms.vstring())
    # # # # # # # # # # for i in scaleSet[1]:
        # # # # # # # # # # getattr(process.lheWeightsTable.weightGroups,scaleSet[0]).append("%i"%(i))


# # # # # # # # # # process.load('RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff')
# # # # # # # # # # process.load('nanotron.NANOProducer.adaptedSV_cff')


# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # Muon Filter

# # # # # # # # # # """
# # # # # # # # # # process.selectedMuonsForFilter = cms.EDFilter("CandViewSelector",
    # # # # # # # # # # src = cms.InputTag("slimmedMuons"),
    # # # # # # # # # # cut = cms.string("pt > 0.0 && isGlobalMuon()")
# # # # # # # # # # )

# # # # # # # # # # process.selectedMuonsMinFilter = cms.EDFilter("CandViewCountFilter",
    # # # # # # # # # # src       = cms.InputTag("selectedMuonsForFilter"),
    # # # # # # # # # # minNumber = cms.uint32(1)
# # # # # # # # # # )

# # # # # # # # # # process.muonFilterSequence = cms.Sequence(
    # # # # # # # # # # process.selectedMuonsForFilter + process.selectedMuonsMinFilter
# # # # # # # # # # )
# # # # # # # # # # """

# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # Electron Filter

# # # # # # # # # # """
# # # # # # # # # # process.selectedElectronsForFilter = cms.EDFilter("CandViewSelector",
    # # # # # # # # # # src = cms.InputTag("slimmedElectrons"),
    # # # # # # # # # # cut = cms.string("pt > 0.0")
# # # # # # # # # # )

# # # # # # # # # # process.selectedElectronsMinFilter = cms.EDFilter("CandViewCountFilter",
    # # # # # # # # # # src       = cms.InputTag("selectedElectronsForFilter"),
    # # # # # # # # # # minNumber = cms.uint32(1)
# # # # # # # # # # )

# # # # # # # # # # process.electronFilterSequence = cms.Sequence(
    # # # # # # # # # # process.selectedElectronsForFilter + process.selectedElectronsMinFilter
# # # # # # # # # # )
# # # # # # # # # # """

# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # Customisation function from PhysicsTools.NanoAOD.nano_cff

# # # # # # # # # # # from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData, nanoAOD_customizeMC
# # # # # # # # # # from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon

# # # # # # # # # # if options.isData:
    # # # # # # # # # # # process = nanoAOD_customizeData(process)
    # # # # # # # # # # process = nanoAOD_customizeCommon(process)
# # # # # # # # # # else:
    # # # # # # # # # # # process = nanoAOD_customizeMC(process)
    # # # # # # # # # # process = nanoAOD_customizeCommon(process)

# # # # # # # # # # # ------------------------------------------------------------------------
# # # # # # # # # # # B-parking muon selection from:
# # # # # # # # # # # https://github.com/DiElectronX/BParkingNANO/blob/main/BParkingNano/python/muonsBPark_cff.py

# # # # # # # # # # Path=["HLT_Mu7_IP4","HLT_Mu8_IP6","HLT_Mu8_IP5","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","HLT_Mu9_IP6","HLT_Mu9_IP5","HLT_Mu9_IP4","HLT_Mu10p5_IP3p5","HLT_Mu12_IP6"]
# # # # # # # # # # #Path=["HLT_Mu9_IP6"]

# # # # # # # # # # process.muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                # # # # # # # # # # muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                # # # # # # # # # # bits = cms.InputTag("TriggerResults","","HLT"),
                                # # # # # # # # # # prescales = cms.InputTag("patTrigger"),
                                # # # # # # # # # # objects = cms.InputTag("slimmedPatTrigger"),
                                # # # # # # # # # # vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),

                                # # # # # # # # # # ##for the output trigger matched collection
                                # # # # # # # # # # maxdR_matching = cms.double(0.1),

                                # # # # # # # # # # ## for the output selected collection (tag + all compatible in dZ)
                                # # # # # # # # # # filterMuon = cms.bool(True),
                                # # # # # # # # # # dzForCleaning_wrtTrgMuon = cms.double(1.0),

                                # # # # # # # # # # ptMin = cms.double(0.5),
                                # # # # # # # # # # absEtaMax = cms.double(2.4),
                                # # # # # # # # # # # keeps only muons with at soft Quality flag
                                # # # # # # # # # # softMuonsOnly = cms.bool(False),
                                # # # # # # # # # # HLTPaths=cms.vstring(Path)#, ### comma to the softMuonsOnly
# # # # # # # # # # #				 L1seeds=cms.vstring(Seed)
                             # # # # # # # # # # )
# # # # # # # # # # #cuts minimun number in B both mu and e, min number of trg, dz muon, dz and dr track, 

# # # # # # # # # # #process.countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
# # # # # # # # # # #    minNumber = cms.uint32(1),
# # # # # # # # # # #    maxNumber = cms.uint32(999999),
# # # # # # # # # # #    src = cms.InputTag("muonTrgSelector", "trgMuons")
# # # # # # # # # # #)


# # # # # # # # # # # process.muonVerticesCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
# # # # # # # # # # #     src = cms.InputTag("muonVerticesTable"),
# # # # # # # # # # #     cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
# # # # # # # # # # #     name = cms.string("muonSV"),
# # # # # # # # # # #     singleton = cms.bool(False), # the number of entries is variable
# # # # # # # # # # #     extension = cms.bool(True), 
# # # # # # # # # # #     variables = cms.PSet(P4Vars,
# # # # # # # # # # #         x   = Var("position().x()", float, doc = "secondary vertex X position, in cm",precision=10),
# # # # # # # # # # #         # y   = Var("position().y()", float, doc = "secondary vertex Y position, in cm",precision=10),
# # # # # # # # # # #         # z   = Var("position().z()", float, doc = "secondary vertex Z position, in cm",precision=14),
# # # # # # # # # # #         # ndof   = Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
# # # # # # # # # # #         # chi2   = Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
# # # # # # # # # # #     ),
# # # # # # # # # # # )
# # # # # # # # # # # process.muonVerticesCandidateTable.variables.pt.precision=10
# # # # # # # # # # # process.muonVerticesCandidateTable.variables.phi.precision=12

# # # # # # # # # # process.muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    # # # # # # # # # # src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    # # # # # # # # # # cut = cms.string(""), #we should not filter on cross linked collections
    # # # # # # # # # # name = cms.string("MuonBPark"),
    # # # # # # # # # # doc  = cms.string("slimmedMuons for BPark after basic selection"),
    # # # # # # # # # # singleton = cms.bool(False), # the number of entries is variable
    # # # # # # # # # # extension = cms.bool(False), # this is the main table for the muons
    # # # # # # # # # # variables = cms.PSet(CandVars,
        # # # # # # # # # # ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        # # # # # # # # # # dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        # # # # # # # # # # dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        # # # # # # # # # # dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        # # # # # # # # # # dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        # # # # # # # # # # vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        # # # # # # # # # # sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
# # # # # # # # # # #        segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
# # # # # # # # # # #        nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
        # # # # # # # # # # #nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
# # # # # # # # # # #        pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        # # # # # # # # # # pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        # # # # # # # # # # pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
# # # # # # # # # # #        tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        # # # # # # # # # # isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
        # # # # # # # # # # isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
        # # # # # # # # # # isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
        # # # # # # # # # # mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
# # # # # # # # # # #        mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        # # # # # # # # # # tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
        # # # # # # # # # # softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
# # # # # # # # # # #        softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
# # # # # # # # # # #        highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        # # # # # # # # # # pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        # # # # # # # # # # tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
# # # # # # # # # # #        mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
# # # # # # # # # # #        miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
# # # # # # # # # # #        multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
        # # # # # # # # # # triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
# # # # # # # # # # #        inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
        # # # # # # # # # # isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco muon is also triggering"),#########################################################3,
# # # # # # # # # # #        toWhichHLTisMatched = Var("userInt('ToWhichHLTisMatched')",int,doc="To which HLT muons is the reco muon matched, -1 for probe" ),
        # # # # # # # # # # matched_dr = Var("userFloat('DR')",float,doc="dr with the matched triggering muon" ),
        # # # # # # # # # # matched_dpt = Var("userFloat('DPT')",float,doc="dpt/pt with the matched triggering muon" ),        #comma
        # # # # # # # # # # skipMuon = Var("userInt('skipMuon')",bool,doc="Is muon skipped (due to large dZ w.r.t. trigger)?"),
        # # # # # # # # # # looseId = Var("userInt('looseId')",int,doc="reco muon is Loose"),
        # # # # # # # # # # fired_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5')",int,doc="reco muon fired this trigger"),
        # # # # # # # # # # fired_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6')",int,doc="reco muon fired this trigger")#,
    # # # # # # # # # # ),
# # # # # # # # # # )

# # # # # # # # # # """
# # # # # # # # # # process.muonsBParkMCMatchForTable = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    # # # # # # # # # # src         = process.muonBParkTable.src,                   # final reco collection
    # # # # # # # # # # matched     = cms.InputTag("finalGenParticlesBPark"),       # final mc-truth particle collection
    # # # # # # # # # # mcPdgId     = cms.vint32(13),                               # one or more PDG ID (13 = mu); absolute values (see below)
    # # # # # # # # # # checkCharge = cms.bool(False),                              # True = require RECO and MC objects to have the same charge
    # # # # # # # # # # mcStatus    = cms.vint32(1),                                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    # # # # # # # # # # maxDeltaR   = cms.double(0.03),                             # Minimum deltaR for the match
    # # # # # # # # # # maxDPtRel   = cms.double(0.5),                              # Minimum deltaPt/Pt for the match
    # # # # # # # # # # resolveAmbiguities    = cms.bool(True),                     # Forbid two RECO objects to match to the same GEN object
    # # # # # # # # # # resolveByMatchQuality = cms.bool(True),                     # False = just match input in order; True = pick lowest deltaR pair first
# # # # # # # # # # )

# # # # # # # # # # process.muonBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    # # # # # # # # # # src     = process.muonBParkTable.src,
    # # # # # # # # # # mcMap   = cms.InputTag("muonsBParkMCMatchForTable"),
    # # # # # # # # # # objName = process.muonBParkTable.name,
    # # # # # # # # # # objType = process.muonBParkTable.name, 
    # # # # # # # # # # branchName = cms.string("genPart"),
    # # # # # # # # # # docString = cms.string("MC matching to status==1 muons"),
# # # # # # # # # # )

# # # # # # # # # # process.selectedMuonsMCMatchEmbedded = cms.EDProducer('MuonMatchEmbedder',
    # # # # # # # # # # src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    # # # # # # # # # # matching = cms.InputTag('muonsBParkMCMatchForTable')
# # # # # # # # # # )
# # # # # # # # # # """

# # # # # # # # # # process.muonTriggerMatchedTable = process.muonBParkTable.clone(
    # # # # # # # # # # src = cms.InputTag("muonTrgSelector:trgMuons"),
    # # # # # # # # # # name = cms.string("TriggerMuon"),
    # # # # # # # # # # doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
    # # # # # # # # # # variables = cms.PSet(CandVars,
        # # # # # # # # # # vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        # # # # # # # # # # ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        # # # # # # # # # # dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        # # # # # # # # # # dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        # # # # # # # # # # dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        # # # # # # # # # # dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        # # # # # # # # # # ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        # # # # # # # # # # sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10)
# # # # # # # # # # #        trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
   # # # # # # # # # # )
# # # # # # # # # # )

# # # # # # # # # # from  PhysicsTools.NanoAOD.triggerObjects_cff import *

# # # # # # # # # # process.triggerObjectBParkTable = cms.EDProducer("TriggerObjectTableBParkProducer",
    # # # # # # # # # # name= cms.string("TrigObjBPark"),
    # # # # # # # # # # src = cms.InputTag("unpackedPatTrigger"),
    # # # # # # # # # # l1Muon = cms.InputTag("gmtStage2Digis","Muon"),
    # # # # # # # # # # selections = cms.VPSet(
        # # # # # # # # # # cms.PSet(
            # # # # # # # # # # name = cms.string("Muon"),
            # # # # # # # # # # id = cms.int32(13),
            # # # # # # # # # # sel = cms.string("type(83) && pt > 5 && coll('hltIterL3MuonCandidates')"), 
            # # # # # # # # # # l1seed = cms.string("type(-81)"), l1deltaR = cms.double(0.5),
            # # # # # # # # # # l2seed = cms.string("type(83) && coll('hltL2MuonCandidates')"),  l2deltaR = cms.double(0.3),
            # # # # # # # # # # qualityBits = cms.string("filter('hltL3fL1s*Park*')"), qualityBitsDoc = cms.string("1 = Muon filters for BPH parking"),
        # # # # # # # # # # ),
    # # # # # # # # # # ),
# # # # # # # # # # )

# # # # # # # # # # # B-parking collection sequences
# # # # # # # # # # process.muonBParkSequence  = cms.Sequence(process.muonTrgSelector)# * process.countTrgMuons)
# # # # # # # # # # process.muonBParkTables    = cms.Sequence(process.muonBParkTable)
# # # # # # # # # # process.muonVertexSequence = cms.Sequence(process.muonVerticesTable)
# # # # # # # # # # process.muonTriggerMatchedTables = cms.Sequence(process.muonTriggerMatchedTable)   ####
# # # # # # # # # # process.triggerObjectBParkTables = cms.Sequence(unpackedPatTrigger + process.triggerObjectBParkTable)
#process.muonBParkMC       = cms.Sequence(process.muonsBParkMCMatchForTable + process.selectedMuonsMCMatchEmbedded + process.muonBParkMCTable)

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
        + process.vertexSequence + process.muonVerticesTable
        # + process.linkedObjects
        # process.TablesTask
        #process.muonFilterSequence+
        # # # # process.nanoSequence+
        # # # # process.adaptedVertexing+
        
        # # # # process.pfOnionTagInfos+
        # # # # process.nanoTable,
        # # # # process.jetTask,
        # process.jetForMETTask 
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
        + process.vertexSequence + process.muonVerticesTable
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
