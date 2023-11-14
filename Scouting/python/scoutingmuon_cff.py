import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer


run3ScoutingMuonRecoTrack = cms.EDProducer("Run3ScoutingMuonRecoTrackProducer",
    muonSource=cms.InputTag("hltScoutingMuonPacker"))

run3ScoutingMuonToPatMuon = cms.EDProducer("Run3ScoutingMuonToPatMuonProducer",
    muonSource=cms.InputTag("hltScoutingMuonPacker"),
    particleSource=cms.InputTag("hltScoutingPFPacker"),
    trackSource=cms.InputTag("run3ScoutingMuonRecoTrack")
)

ScoutingMuonTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
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
            expr = cms.string("userInt('numberofmatchedstations')"),
            #expr = cms.string("1"),
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
        # scouting variables
        isPFmatched = cms.PSet(
            doc = cms.string('muon is PF matched'),
            expr = cms.string("userInt('isPFmatched')"),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        normalizedChi2 = cms.PSet(
            doc = cms.string('normalizedChi2'),
            expr = cms.string("userFloat('normalizedChi2')"),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        ecalIso = cms.PSet(
            doc = cms.string('ecalIso'),
            expr = cms.string("userFloat('ecalIso')"),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        hcalIso = cms.PSet(
            doc = cms.string('hcalIso'),
            expr = cms.string("userFloat('hcalIso')"),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        trackIso = cms.PSet(
            doc = cms.string('trackIso'),
            expr = cms.string("userFloat('trackIso')"),
            precision = cms.int32(20),
            type = cms.string('float')
        ),
        nValidStandAloneMuonHits = cms.PSet(
            doc = cms.string('nValidStandAloneMuonHits'),
            expr = cms.string("userInt('nValidStandAloneMuonHits')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nStandAloneMuonMatchedStations = cms.PSet(
            doc = cms.string('nStandAloneMuonMatchedStations'),
            expr = cms.string("userInt('nStandAloneMuonMatchedStations')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nValidRecoMuonHits = cms.PSet(
            doc = cms.string('nValidRecoMuonHits'),
            expr = cms.string("userInt('nValidRecoMuonHits')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonChambers = cms.PSet(
            doc = cms.string('nRecoMuonChambers'),
            expr = cms.string("userInt('nRecoMuonChambers')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonChambersCSCorDT = cms.PSet(
            doc = cms.string('nRecoMuonChambersCSCorDT'),
            expr = cms.string("userInt('nRecoMuonChambersCSCorDT')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonMatches = cms.PSet(
            doc = cms.string('nRecoMuonMatches'),
            expr = cms.string("userInt('nRecoMuonMatches')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonMatchedStations = cms.PSet(
            doc = cms.string('nRecoMuonMatchedStations'),
            expr = cms.string("userInt('nRecoMuonMatchedStations')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonExpectedMatchedStations = cms.PSet(
            doc = cms.string('nRecoMuonExpectedMatchedStations'),
            expr = cms.string("userInt('nRecoMuonExpectedMatchedStations')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        recoMuonStationMask = cms.PSet(
            doc = cms.string('recoMuonStationMask'),
            expr = cms.string("userInt('recoMuonStationMask')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nRecoMuonMatchedRPCLayers = cms.PSet(
            doc = cms.string('nRecoMuonMatchedRPCLayers'),
            expr = cms.string("userInt('nRecoMuonMatchedRPCLayers')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        recoMuonRPClayerMask = cms.PSet(
            doc = cms.string('recoMuonRPClayerMask'),
            expr = cms.string("userInt('recoMuonRPClayerMask')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nValidPixelHits = cms.PSet(
            doc = cms.string('nValidPixelHits'),
            expr = cms.string("userInt('nValidPixelHits')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nValidStripHits = cms.PSet(
            doc = cms.string('nValidStripHits'),
            expr = cms.string("userInt('nValidStripHits')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nPixelLayersWithMeasurement = cms.PSet(
            doc = cms.string('nPixelLayersWithMeasurement'),
            expr = cms.string("userInt('nPixelLayersWithMeasurement')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nTrackerLayersWithMeasurement = cms.PSet(
            doc = cms.string('nTrackerLayersWithMeasurement'),
            expr = cms.string("userInt('nTrackerLayersWithMeasurement')"),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
    )
)

#This sequence produces PAT Muons out of Scouting muons
muonSequence = cms.Sequence(run3ScoutingMuonRecoTrack * run3ScoutingMuonToPatMuon * ScoutingMuonTable)