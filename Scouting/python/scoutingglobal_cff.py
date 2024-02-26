#This adds global variables like rho and MET
import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer


#These are imported from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_run3scouting_cff.py (by A. Lintuluoto)
rhoScoutingTable = cms.EDProducer("GlobalVariablesTableProducer",
    name = cms.string(""),
    variables = cms.PSet(
        ScoutingRho = ExtVar( cms.InputTag("hltScoutingPFPacker", "rho"), "double", doc = "rho from all scouting PF Candidates, used e.g. for JECs" ),
    )
)

metScoutingTable = cms.EDProducer("GlobalVariablesTableProducer",
    name = cms.string("MET"),
    variables = cms.PSet(
        pt = ExtVar( cms.InputTag("hltScoutingPFPacker", "pfMetPt"), "double", doc = "scouting MET pt"),
        phi = ExtVar( cms.InputTag("hltScoutingPFPacker", "pfMetPhi"), "double", doc = "scouting MET phi"),
    )
)

globalSequence = cms.Sequence(rhoScoutingTable * metScoutingTable)