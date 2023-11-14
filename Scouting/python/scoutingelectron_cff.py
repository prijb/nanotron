import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

run3ScoutingEleToPatEle = cms.EDProducer("Run3ScoutingEleToPatEleProducer",
    eleSource=cms.InputTag("hltScoutingEgammaPacker"))

run3ScoutingPhotonToPatPhoton = cms.EDProducer("Run3ScoutingPhotonToPatPhotonProducer",
    photonSource=cms.InputTag("hltScoutingEgammaPacker"))

ScoutingElectronTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
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

electronSequence = cms.Sequence((run3ScoutingEleToPatEle * ScoutingElectronTable) + run3ScoutingPhotonToPatPhoton)

