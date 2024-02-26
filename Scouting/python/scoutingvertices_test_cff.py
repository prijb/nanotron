import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

#This provides the PVs that are used for the inclusive muon vertexing
run3ScoutingVertices = cms.EDProducer("Run3ScoutingVtxToVtxProducer",
    pvSource=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    svSource=cms.InputTag("hltScoutingMuonPacker", "displacedVtx"),
    muonSource=cms.InputTag("hltScoutingMuonPacker")
)

pvScoutingTable = cms.EDProducer("SimpleRun3ScoutingVertexFlatTableProducer",
     src = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
     cut = cms.string(""),
     name = cms.string("PV"),
     doc  = cms.string("PrimaryVertex scouting information"),
     singleton = cms.bool(False),
     extension = cms.bool(False),
     variables = cms.PSet(
         x = Var('x', 'float', precision=10, doc='position x coordinate'),
         y = Var('y', 'float', precision=10, doc='position y coordinate'),
         z = Var('z', 'float', precision=10, doc='position z coordinate'),
         xError = Var('xError', 'float', precision=10, doc='x error'),
         yError = Var('yError', 'float', precision=10, doc='y error'),
         zError = Var('zError', 'float', precision=10, doc='z error'),
         tracksSize = Var('tracksSize', 'int', doc='number of tracks'),
         chi2 = Var('chi2', 'float', precision=10, doc='chi squared'),
         ndof = Var('ndof', 'int', doc='number of degrees of freedom'),
         isValidVtx = Var('isValidVtx', 'bool', doc='is valid'),
     )
)

svScoutingTable = cms.EDProducer("Run3ScoutingVtxTableProducer",
    pvSource=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    svSource=cms.InputTag("hltScoutingMuonPacker", "displacedVtx"),
    muonSource=cms.InputTag("hltScoutingMuonPacker"),
    svName  = cms.string("SV")
)

vertexSequence = cms.Sequence(run3ScoutingVertices * pvScoutingTable * svScoutingTable)