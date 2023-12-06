import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

run3ScoutingJetToPatJet = cms.EDProducer("Run3ScoutingJetToPatJetProducer",
    jetSource=cms.InputTag("hltScoutingPFPacker"))

ScoutingJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string('slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15)'),
    extension = cms.bool(False),
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

jetSequence = cms.Sequence(run3ScoutingJetToPatJet * ScoutingJetTable)
