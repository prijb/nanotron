import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.run3scouting_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

#These are imported from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_run3scouting_cff.py (by A. Lintuluoto)
particleTask = cms.Task(scoutingPFCands)

#ak4JetTableTask = cms.Task(ak4ScoutingJets,ak4ScoutingJetParticleNetJetTagInfos,ak4ScoutingJetParticleNetJetTags,ak4ScoutingJetTable)
#ak8JetTableTask = cms.Task(ak8ScoutingJets,ak8ScoutingJetsSoftDrop,ak8ScoutingJetsSoftDropMass,ak8ScoutingJetEcfNbeta1,ak8ScoutingJetNjettiness,ak8ScoutingJetParticleNetJetTagInfos,ak8ScoutingJetParticleNetJetTags,ak8ScoutingJetParticleNetMassRegressionJetTags,ak8ScoutingJetTable)

ak4RenamedScoutingJetTable = ak4ScoutingJetTable.clone(
    name = cms.string("Jet"),
    doc = cms.string("ScoutingJet"),
)

ak8RenamedScoutingJetTable = ak8ScoutingJetTable.clone(
    name = cms.string("FatJet"),
    doc = cms.string("ScoutingFatJet"),
)

ak4JetTableTask = cms.Task(ak4ScoutingJets,ak4ScoutingJetParticleNetJetTagInfos,ak4ScoutingJetParticleNetJetTags,ak4RenamedScoutingJetTable)
ak8JetTableTask = cms.Task(ak8ScoutingJets,ak8ScoutingJetsSoftDrop,ak8ScoutingJetsSoftDropMass,ak8ScoutingJetEcfNbeta1,ak8ScoutingJetNjettiness,ak8ScoutingJetParticleNetJetTagInfos,ak8ScoutingJetParticleNetJetTags,ak8ScoutingJetParticleNetMassRegressionJetTags,ak8RenamedScoutingJetTable)

jetSequence = cms.Sequence(particleTask, ak4JetTableTask, ak8JetTableTask)