import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

run3ScoutingVertices = cms.EDProducer("Run3ScoutingVtxToVtxProducer",
    pvSource=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    svSource=cms.InputTag("hltScoutingMuonPacker", "displacedVtx"),
    muonSource=cms.InputTag("hltScoutingMuonPacker")
)

ScoutingVertexTable = cms.EDProducer("VertexTableProducer",
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

ScoutingSVTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
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
    src = cms.InputTag("ScoutingVertexTable"),
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

vertexSequence = cms.Sequence(run3ScoutingVertices * ScoutingVertexTable * (ScoutingSVTable))