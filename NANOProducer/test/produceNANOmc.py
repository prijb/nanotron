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
process.MessageLogger.cerr.FwkReport.reportEvery = 100
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
               ['llpnanoAOD_step_mu'] if options.isData else ['llpnanoAOD_step'] # ['llpnanoAOD_step_mu','llpnanoAOD_step_ele'] ~ boolean OR (union) between 'mu' and 'ele' paths
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



process.load('nanotron.NANOProducer.GenDisplacedVertices_cff')

# ------------------------------------------------------------------------
# Parton Density Functions

process.load('RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff')
process.load('nanotron.NANOProducer.adaptedSV_cff')

# ------------------------------------------------------------------------
# B-parking muon selection from:
# https://github.com/DiElectronX/BParkingNANO/blob/main/BParkingNano/python/muonsBPark_cff.py

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

#Path=["HLT_Mu9_IP6"]

process.muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                bits = cms.InputTag("TriggerResults","","HLT"),
                                prescales = cms.InputTag("patTrigger"),
                                objects = cms.InputTag("slimmedPatTrigger"),
                                vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),

                                ##for the output trigger matched collection
                                maxdR_matching = cms.double(0.1),

                                ## for the output selected collection (tag + all compatible in dZ)
                                filterMuon = cms.bool(True),
                                dzForCleaning_wrtTrgMuon = cms.double(1.0),

                                ptMin = cms.double(0.5),
                                absEtaMax = cms.double(2.4),
                                # keeps only muons with at soft Quality flag
                                softMuonsOnly = cms.bool(False),
                                HLTPaths=cms.vstring(Path)#, ### comma to the softMuonsOnly
#				 L1seeds=cms.vstring(Seed)
                             )
#cuts minimun number in B both mu and e, min number of trg, dz muon, dz and dr track, 

#process.countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
#    minNumber = cms.uint32(1),
#    maxNumber = cms.uint32(999999),
#    src = cms.InputTag("muonTrgSelector", "trgMuons")
#)

process.muonVerticesTable = cms.EDProducer("MuonVertexProducer",
    srcMuon = cms.InputTag("slimmedMuons"),
    pvSrc   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    ptMin   = cms.double(0.8),
    svName  = cms.string("muonSV"),
)

process.fourmuonVerticesTable = cms.EDProducer("FourMuonVertexProducer",
    srcMuon = cms.InputTag("slimmedMuons"),
    #srcMuon = cms.InputTag("finalMuons"),
    pvSrc   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    svCut   = cms.string(""),  # careful: adding a cut here would make the collection matching inconsistent with the SV table
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    ptMin   = cms.double(0.8),
    svName  = cms.string("fourmuonSV"),
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
        skipMuon = Var("userInt('skipMuon')",bool,doc="Is muon skipped (due to large dZ w.r.t. trigger)?"),
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
    ),
)

# B-parking collection sequences
process.muonBParkSequence  = cms.Sequence(process.muonTrgSelector)# * process.countTrgMuons)
process.muonBParkTables    = cms.Sequence(process.muonBParkTable)
#process.muonVertexSequence = cms.Sequence(process.muonVerticesTable)
process.muonVertexSequence = cms.Sequence(process.muonVerticesTable + process.fourmuonVerticesTable)
process.muonTriggerMatchedTables = cms.Sequence(process.muonTriggerMatchedTable)   ####
process.triggerObjectBParkTables = cms.Sequence(unpackedPatTrigger + process.triggerObjectBParkTable + process.triggerMuonTable)
#process.muonBParkMC       = cms.Sequence(process.muonsBParkMCMatchForTable + process.selectedMuonsMCMatchEmbedded + process.muonBParkMCTable)

# ------------------------------------------------------------------------

# ========================================================================
# ** DATA SEQUENCE **
# ========================================================================

if options.isData:

    # Main
    process.llpnanoAOD_step_mu = cms.Path(
        process.nanoSequence+
        process.adaptedVertexing
    )



    # B-parking additions
    process.llpnanoAOD_step_mu += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence
    #process.llpnanoAOD_step_ele += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence
    #process.llpnanoAOD_step_mu += process.metadata

# ========================================================================
# ** MC SEQUENCE **
# ========================================================================

else:
    # Main
    process.llpnanoAOD_step = cms.Path(
        process.nanoSequenceMC
    )

    # B-parking additions
    process.llpnanoAOD_step += process.muonBParkSequence + process.muonBParkTables + process.muonTriggerMatchedTables + process.triggerObjectBParkTables + process.muonVertexSequence

process.endjob_step           = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# ------------------------------------------------------------------------
# Sequence

if options.isData:
#    process.schedule = cms.Schedule(process.llpnanoAOD_step_mu, process.llpnanoAOD_step_ele, process.endjob_step, process.NANOAODSIMoutput_step)
    process.schedule = cms.Schedule(process.llpnanoAOD_step_mu, process.endjob_step, process.NANOAODSIMoutput_step)
else:
    process.schedule = cms.Schedule(process.llpnanoAOD_step, process.endjob_step, process.NANOAODSIMoutput_step)

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
    
    "l1bits",
]

#override final jets

#process.finalJets.addBTagInfo=cms.bool(True)
#process.finalJets.addDiscriminators = cms.bool(True)
#process.finalJets.addTagInfos=cms.bool(True)

for moduleName in modulesToRemove:
    if hasattr(process,moduleName):
        print("removing module:", moduleName)
        if options.isData:
            process.nanoSequence.remove(getattr(process,moduleName))
        else:
            process.nanoSequenceMC.remove(getattr(process,moduleName))
    else:
        print("module for removal not found:", moduleName)

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