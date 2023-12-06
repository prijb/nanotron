import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer


run3ScoutingParticles = cms.EDProducer("Run3ScoutingPFToCandidateProducer",
    vertexSource=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    particleSource=cms.InputTag("hltScoutingPFPacker")
)

ScoutingPfparticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string(''),
    extension = cms.bool(False),
    externalVariables = cms.PSet(

    ),
    maxLen = cms.optional.uint32,
    mightGet = cms.optional.untracked.vstring,
    name = cms.string('pfpart'),
    singleton = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    # src = cms.InputTag("run3ScoutingVertices", "svs"),
    src = cms.InputTag("run3ScoutingParticles"),
    variables = cms.PSet(
        eta = cms.PSet(
            doc = cms.string('eta'),
            expr = cms.string('eta'),
            precision = cms.int32(12),
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
        x = cms.PSet(
            doc = cms.string('primary vertex X position, in cm'),
            expr = cms.string('vx'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        y = cms.PSet(
            doc = cms.string('primary vertex Y position, in cm'),
            expr = cms.string('vy'),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        z = cms.PSet(
            doc = cms.string('primary vertex Z position, in cm'),
            expr = cms.string('vz'),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        pdgId = cms.PSet(
            doc = cms.string('pdgId'),
            expr = cms.string('pdgId'),
            precision = cms.int32(-1),
            type = cms.string('int')
        )
    )
)

particleSequence = cms.Sequence(run3ScoutingParticles * ScoutingPfparticleTable)