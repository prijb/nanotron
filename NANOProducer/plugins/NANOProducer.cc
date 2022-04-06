//
//
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "nanotron/DataFormats/interface/OnionTagInfo.h"
#include "nanotron/DataFormats/interface/MCLabel.h"
#include "nanotron/DataFormats/interface/MCLabelInfo.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <Math/Vector4D.h>

#include "FlatTableFiller.h"

//
// class declaration
//

class NANOProducer : public edm::stream::EDProducer<> {
    public:
        explicit NANOProducer(const edm::ParameterSet&);
        ~NANOProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        PropertyList<nanotron::JetFeatures> globalProperties;
        PropertyList<nanotron::ShallowTagInfoFeatures> csvProperties;
        PropertyList<nanotron::ChargedCandidateFeatures> cpfProperties;
        PropertyList<nanotron::NeutralCandidateFeatures> npfProperties;
        PropertyList<nanotron::SecondaryVertexFeatures> svProperties;
        PropertyList<nanotron::MuonCandidateFeatures> muonProperties;
        PropertyList<nanotron::ElectronCandidateFeatures> electronProperties;
        
   private:
      const edm::EDGetTokenT<edm::View<pat::Jet>> _jet_src;
      const edm::EDGetTokenT<std::vector<reco::OnionTagInfo>> _tag_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOProducer::NANOProducer(const edm::ParameterSet& iConfig) :
    _jet_src(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    _tag_src(consumes<std::vector<reco::OnionTagInfo>>(iConfig.getParameter<edm::InputTag>("srcTags"))) {
    
    globalProperties = {

        PROPERTY(nanotron::JetFeatures, pt,      "doc"),
        PROPERTY(nanotron::JetFeatures, eta,     "doc"),
        PROPERTY(nanotron::JetFeatures, phi,     "doc"),
        PROPERTY(nanotron::JetFeatures, mass,    "doc"),
        PROPERTY(nanotron::JetFeatures, energy,  "doc"),

        PROPERTY(nanotron::JetFeatures, area,    "doc"),

        PROPERTY(nanotron::JetFeatures, beta,    "doc"),
        PROPERTY(nanotron::JetFeatures, dR2Mean, "doc"),
        PROPERTY(nanotron::JetFeatures, frac01,  "doc"),
        PROPERTY(nanotron::JetFeatures, frac02,  "doc"),
        PROPERTY(nanotron::JetFeatures, frac03,  "doc"),
        PROPERTY(nanotron::JetFeatures, frac04,  "doc"),

        PROPERTY(nanotron::JetFeatures, jetR,    "doc"),
        PROPERTY(nanotron::JetFeatures, jetRchg, "doc"),

        PROPERTY(nanotron::JetFeatures, n60, "doc"),
        PROPERTY(nanotron::JetFeatures, n90, "doc"),

        PROPERTY(nanotron::JetFeatures, chargedEmEnergyFraction,     "doc"),
        PROPERTY(nanotron::JetFeatures, chargedHadronEnergyFraction, "doc"),
        PROPERTY(nanotron::JetFeatures, chargedMuEnergyFraction,     "doc"),
        PROPERTY(nanotron::JetFeatures, electronEnergyFraction,      "doc"),

        PROPERTY(nanotron::JetFeatures, tau1, "doc"),
        PROPERTY(nanotron::JetFeatures, tau2, "doc"),
        PROPERTY(nanotron::JetFeatures, tau3, "doc"),

        PROPERTY(nanotron::JetFeatures, relMassDropMassAK, "doc"),
        PROPERTY(nanotron::JetFeatures, relMassDropMassCA, "doc"),
        PROPERTY(nanotron::JetFeatures, relSoftDropMassAK, "doc"),
        PROPERTY(nanotron::JetFeatures, relSoftDropMassCA, "doc"),

        PROPERTY(nanotron::JetFeatures, thrust,          "doc"),
        PROPERTY(nanotron::JetFeatures, sphericity,      "doc"),
        PROPERTY(nanotron::JetFeatures, circularity,     "doc"),
        PROPERTY(nanotron::JetFeatures, isotropy,        "doc"),
        PROPERTY(nanotron::JetFeatures, eventShapeC,     "doc"),
        PROPERTY(nanotron::JetFeatures, eventShapeD,     "doc"),

        PROPERTY(nanotron::JetFeatures, numberCpf,       "doc"),
        PROPERTY(nanotron::JetFeatures, numberNpf,       "doc"),
        PROPERTY(nanotron::JetFeatures, numberSv,        "doc"),
        PROPERTY(nanotron::JetFeatures, numberSvAdapted, "doc"),
        PROPERTY(nanotron::JetFeatures, numberMuon,      "doc"),
        PROPERTY(nanotron::JetFeatures, numberElectron,  "doc")
    };
    
    csvProperties = {
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSumJetEtRatio,      "ratio of track sum transverse energy over jet energy"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSumJetDeltaR,       "pseudoangular distance between jet axis and track fourvector sum"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, vertexCategory,          "category of secondary vertex (Reco, Pseudo, No)"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSip2dValAboveCharm, "track 2D signed impact parameter of first track lifting mass above charm"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSip2dSigAboveCharm, "track 2D signed impact parameter significance of first track lifting mass above charm"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSip3dValAboveCharm, "track 3D signed impact parameter of first track lifting mass above charm"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, trackSip3dSigAboveCharm, "track 3D signed impact parameter significance of first track lifting mass above charm"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, jetNTracksEtaRel,        "tracks associated to jet for which trackEtaRel is calculated"),
        PROPERTY(nanotron::ShallowTagInfoFeatures, jetNSelectedTracks,      "doc")
    };
    
    cpfProperties = {
        PROPERTY(nanotron::ChargedCandidateFeatures, ptrel,  "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, deta,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, dphi,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, deltaR, "doc"),

        PROPERTY(nanotron::ChargedCandidateFeatures, px, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, py, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, pz, "doc"),

        PROPERTY(nanotron::ChargedCandidateFeatures, trackEtaRel, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackPtRel,  "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackPPar,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackDeltaR, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackPParRatio,  "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackPtRatio,    "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dVal,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dSig,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dVal,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dSig,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackJetDistVal, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, trackJetDistSig, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, drminsv, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, vertex_association, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, fromPV,  "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, puppi_weight,  "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_chi2,    "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_quality, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_numberOfValidPixelHits,     "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_pixelLayersWithMeasurement, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_numberOfValidStripHits,     "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, track_stripLayersWithMeasurement, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, relmassdrop, "doc"),
        
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dValSV, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dSigSV, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dValSV, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dSigSV, "doc"), 

        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dValSV_adapted, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip2dSigSV_adapted, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dValSV_adapted, "doc"), 
        PROPERTY(nanotron::ChargedCandidateFeatures, trackSip3dSigSV_adapted, "doc"), 

        PROPERTY(nanotron::ChargedCandidateFeatures, matchedMuon,       "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, matchedElectron,   "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, matchedSV,         "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, matchedSV_adapted, "doc"),
        
        PROPERTY(nanotron::ChargedCandidateFeatures, track_ndof, "doc"),
        PROPERTY(nanotron::ChargedCandidateFeatures, dZmin, "doc")
    };


    npfProperties = {
        PROPERTY(nanotron::NeutralCandidateFeatures, ptrel,  "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, deta,   "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, dphi,   "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, deltaR, "doc"),

        PROPERTY(nanotron::NeutralCandidateFeatures, px, "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, py, "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, pz, "doc"),

        PROPERTY(nanotron::NeutralCandidateFeatures, isGamma,       "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, hcal_fraction, "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, drminsv,       "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, puppi_weight,  "doc"),
        PROPERTY(nanotron::NeutralCandidateFeatures, relmassdrop,   "doc")
    };
    
    svProperties = {
        PROPERTY(nanotron::SecondaryVertexFeatures, ptrel,   "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, deta,    "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, dphi,    "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, deltaR,  "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, mass,    "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, ntracks, "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, chi2,    "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, ndof,    "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, dxy,     "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, dxysig,  "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, d3d,     "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, d3dsig,  "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, costhetasvpv, "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, enratio,      "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, vx, "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, vy, "doc"),
        PROPERTY(nanotron::SecondaryVertexFeatures, vz, "doc")
    };
    
    muonProperties = {
        PROPERTY(nanotron::MuonCandidateFeatures, isGlobal,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, isTight,      "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, isMedium,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, isLoose,      "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, isStandAlone, "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, ptrel,  "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, deta,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dphi,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, px,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, py,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, pz,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, charge, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, energy, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, et,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, deltaR, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, numberOfMatchedStations, "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, IP2d,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, IP2dSig,  "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, IP3d,     "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, IP3dSig,  "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, EtaRel,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dxy,      "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dxyError, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dxySig,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dz,       "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dzError,  "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, dzSig,    "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, numberOfValidPixelHits, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, numberOfpixelLayersWithMeasurement, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, numberOfstripLayersWithMeasurement, "doc"), //that does not help. needs to be discussed.

        PROPERTY(nanotron::MuonCandidateFeatures, chi2, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, ndof, "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, caloIso, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, ecalIso, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, hcalIso, "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, sumPfChHadronPt,    "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, sumPfNeuHadronEt,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, Pfpileup,           "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, sumPfPhotonEt,      "doc"),

        PROPERTY(nanotron::MuonCandidateFeatures, sumPfChHadronPt03,  "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, sumPfNeuHadronEt03, "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, Pfpileup03,         "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, sumPfPhotonEt03,    "doc"),
        
        PROPERTY(nanotron::MuonCandidateFeatures, timeAtIpInOut,      "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, timeAtIpInOutErr,   "doc"),
        PROPERTY(nanotron::MuonCandidateFeatures, timeAtIpOutIn,      "doc")
    };
    
    
    electronProperties = {
        PROPERTY(nanotron::ElectronCandidateFeatures, ptrel,  "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaR, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deta,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dphi,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, px,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, py,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, pz,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, charge, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, energy, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, EtFromCaloEn, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, isEB,   "doc"), 
        PROPERTY(nanotron::ElectronCandidateFeatures, isEE,   "doc"), 
        PROPERTY(nanotron::ElectronCandidateFeatures, ecalEnergy, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, isPassConversionVeto, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, convDist,       "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, convFlags,      "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, convRadius,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, hadronicOverEm, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, ecalDrivenSeed, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, IP2d,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, IP2dSig, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, IP3d,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, IP3dSig, "doc"),

        PROPERTY(nanotron::ElectronCandidateFeatures, elecSC_energy, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, elecSC_deta,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, elecSC_dphi,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, elecSC_et,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, elecSC_eSuperClusterOverP, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, scPixCharge,       "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, superClusterFbrem, "doc"),

        PROPERTY(nanotron::ElectronCandidateFeatures, eSeedClusterOverP,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, eSeedClusterOverPout, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, eSuperClusterOverP,   "doc"),

        // shower shape
        PROPERTY(nanotron::ElectronCandidateFeatures, sigmaEtaEta,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, sigmaIetaIeta,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, sigmaIphiIphi,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, e5x5,            "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, e5x5Rel,         "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, e1x5Overe5x5,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, e2x5MaxOvere5x5,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, r9,                 "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, hcalOverEcal,       "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, hcalDepth1OverEcal, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, hcalDepth2OverEcal, "doc"),

        // Track-Cluster Matching Attributes
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaEtaEleClusterTrackAtCalo,  "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaEtaSeedClusterTrackAtCalo, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaPhiSeedClusterTrackAtCalo, "doc"), 
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaEtaSeedClusterTrackAtVtx,  "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaEtaSuperClusterTrackAtVtx, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaPhiEleClusterTrackAtCalo,  "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, deltaPhiSuperClusterTrackAtVtx, "doc"),

        PROPERTY(nanotron::ElectronCandidateFeatures, sCseedEta, "doc"),

        // electron gsf variables. 
        PROPERTY(nanotron::ElectronCandidateFeatures, EtaRel,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dxy,      "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dxyError, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dxySig,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dz,       "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dzError,  "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dzSig,    "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, nbOfMissingHits, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, gsfCharge,       "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, ndof,            "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, chi2,            "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, numberOfBrems,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, fbrem,           "doc"),

        // Isolation block
        PROPERTY(nanotron::ElectronCandidateFeatures, neutralHadronIso,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, particleIso,          "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, photonIso,            "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, puChargedHadronIso,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, trackIso,             "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, ecalPFClusterIso,     "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, hcalPFClusterIso,     "doc"),

        PROPERTY(nanotron::ElectronCandidateFeatures, pfSumPhotonEt,        "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, pfSumChargedHadronPt, "doc"), 
        PROPERTY(nanotron::ElectronCandidateFeatures, pfSumNeutralHadronEt, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, pfSumPUPt,            "doc"),

        PROPERTY(nanotron::ElectronCandidateFeatures, dr04TkSumPt,          "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04EcalRecHitSumEt,        "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalDepth1TowerSumEt,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalDepth1TowerSumEtBc, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalDepth2TowerSumEt,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalDepth2TowerSumEtBc, "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalTowerSumEt,   "doc"),
        PROPERTY(nanotron::ElectronCandidateFeatures, dr04HcalTowerSumEtBc, "doc")
    };

    produces<nanoaod::FlatTable>("global");
    produces<nanoaod::FlatTable>("csv");
    produces<nanoaod::FlatTable>("cpf");
    produces<nanoaod::FlatTable>("npf");
    produces<nanoaod::FlatTable>("sv");
    produces<nanoaod::FlatTable>("svAdapted");
    produces<nanoaod::FlatTable>("length");
    produces<nanoaod::FlatTable>("muon");
    produces<nanoaod::FlatTable>("electron");
}


NANOProducer::~NANOProducer() {
}


// ------------ method called to produce the data  ------------
//
void
NANOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(_jet_src, jets);
    edm::Handle<std::vector<reco::OnionTagInfo>> tag_infos;
    iEvent.getByToken(_tag_src, tag_infos);

    unsigned int ntags = tag_infos->size();

    std::vector<int> global_jetIdx;
    std::vector<int> csv_jetIdx;
    std::vector<int> cpf_jetIdx;
    std::vector<int> npf_jetIdx;
    std::vector<int> sv_jetIdx;
    std::vector<int> sv_adapted_jetIdx;
    std::vector<int> mu_jetIdx;
    std::vector<int> elec_jetIdx;

    auto lengthTable = std::make_unique<nanoaod::FlatTable>(ntags, "length", false, false);
    std::vector<int> cpf_length;
    std::vector<int> npf_length;
    std::vector<int> sv_length;
    std::vector<int> sv_adapted_length;

    std::vector<int> elec_length;
    std::vector<int> mu_length;

    auto globalTable = std::make_unique<nanoaod::FlatTable>(ntags, "global", false, false);
    auto csvTable = std::make_unique<nanoaod::FlatTable>(ntags, "csv", false, false);


    FlatTableFillerList<nanotron::JetFeatures> globalFillerList(globalProperties);
    FlatTableFillerList<nanotron::ShallowTagInfoFeatures> csvFillerList(csvProperties);
    
    FlatTableFillerList<nanotron::ChargedCandidateFeatures> cpfFillerList(cpfProperties);
    FlatTableFillerList<nanotron::NeutralCandidateFeatures> npfFillerList(npfProperties);
    FlatTableFillerList<nanotron::SecondaryVertexFeatures> svFillerList(svProperties);
    FlatTableFillerList<nanotron::SecondaryVertexFeatures> svAdaptedFillerList(svProperties);
    FlatTableFillerList<nanotron::MuonCandidateFeatures> muonFillerList(muonProperties);
    FlatTableFillerList<nanotron::ElectronCandidateFeatures> electronFillerList(electronProperties);

    unsigned int nmu_total   = 0;
    unsigned int nelec_total = 0;
    unsigned int ncpf_total  = 0;
    unsigned int nnpf_total  = 0;
    unsigned int nsv_total   = 0;
    unsigned int nsv_total_adapted = 0;

    for (unsigned int itag = 0; itag < ntags; ++itag) {

        const auto& features = tag_infos->at(itag).features();

        unsigned int nmu   = features.mu_features.size();
        unsigned int nelec = features.elec_features.size();
        unsigned int ncpf  = features.cpf_features.size();
        unsigned int nnpf  = features.npf_features.size();
        unsigned int nsv   = features.sv_features.size();
        unsigned int nsv_adapted = features.sv_adapted_features.size();
        
        nmu_total   += nmu;
        nelec_total += nelec;
        ncpf_total  += ncpf;
        nnpf_total  += nnpf;
        nsv_total   += nsv;
        nsv_total_adapted += nsv_adapted;

        cpf_length.push_back(ncpf);
        npf_length.push_back(nnpf);
        sv_length.push_back(nsv);
        sv_adapted_length.push_back(nsv_adapted);
        elec_length.push_back(nelec);
        mu_length.push_back(nmu);

        int jetIdx = -1;
        auto base_jet_ref = tag_infos->at(itag).jet();

        if (base_jet_ref.isAvailable() and base_jet_ref.isNonnull()) {
            const auto& base_jet = base_jet_ref.get();

            for (std::size_t ijet = 0; ijet < jets->size(); ++ijet) {
                auto jet = jets->at(ijet);
                if (reco::deltaR(base_jet->p4(), jet.p4()) < 1e-4) {
                    jetIdx = ijet;
                    break;
                }
            }
        }
        
        global_jetIdx.push_back(jetIdx);
        csv_jetIdx.push_back(jetIdx);
        
        globalFillerList.push_back(features.jet_features);
        csvFillerList.push_back(features.tag_info_features);
    }

    auto muonTable       = std::make_unique<nanoaod::FlatTable>(nmu_total,   "muon", false, false);
    auto electronTable   = std::make_unique<nanoaod::FlatTable>(nelec_total, "electron", false, false);
    auto cpfTable        = std::make_unique<nanoaod::FlatTable>(ncpf_total,  "cpf", false, false);
    auto npfTable        = std::make_unique<nanoaod::FlatTable>(nnpf_total,  "npf", false, false);
    auto svTable         = std::make_unique<nanoaod::FlatTable>(nsv_total,   "sv",  false, false);
    auto svTable_adapted = std::make_unique<nanoaod::FlatTable>(nsv_total_adapted, "svAdapted", false, false);

    for (unsigned int itag = 0; itag < ntags; ++itag) {

        const auto& features = tag_infos->at(itag).features();
        auto mu   = features.mu_features;
        auto elec = features.elec_features;
        auto cpf  = features.cpf_features;
        auto npf  = features.npf_features;
        auto sv   = features.sv_features;
        auto sv_adapted = features.sv_adapted_features;

        unsigned int nmu   = features.mu_features.size();
        unsigned int nelec = features.elec_features.size();
        unsigned int ncpf  = features.cpf_features.size();
        unsigned int nnpf  = features.npf_features.size();
        unsigned int nsv   = features.sv_features.size();
        unsigned int nsv_adapted = features.sv_adapted_features.size();

        int jetIdx = global_jetIdx[itag];

        for (unsigned int i = 0; i < ncpf; ++i) {

            cpf_jetIdx.push_back(jetIdx);
            cpfFillerList.push_back(cpf.at(i));
        }

        for (unsigned int i = 0; i < nnpf; ++i) {

            npf_jetIdx.push_back(jetIdx);          
            npfFillerList.push_back(npf.at(i));
        }

        for (unsigned int i = 0; i < nsv; ++i) {

            sv_jetIdx.push_back(jetIdx);
            svFillerList.push_back(sv.at(i));
        }

        for (unsigned int i = 0; i < nsv_adapted; ++i) {

            sv_adapted_jetIdx.push_back(jetIdx);
            svAdaptedFillerList.push_back(sv_adapted.at(i));
        }

        for (unsigned int i = 0; i < nmu; ++i) {

            mu_jetIdx.push_back(jetIdx);
            muonFillerList.push_back(mu.at(i));
        }

        for (unsigned int i = 0; i < nelec; ++i) {

            elec_jetIdx.push_back(jetIdx);
            electronFillerList.push_back(elec.at(i));
        }
    }

    lengthTable->addColumn<int>("cpf",       cpf_length,        "charged PF candidate track offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("npf",       npf_length,        "neutral PF candidate offset",  nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("sv",        sv_length,         "secondary vertex (SV) offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("svAdapted", sv_adapted_length, "secondary vertex (SV) offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("mu",        mu_length,         "muon offset",      nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("ele",       elec_length,       "electron offset",  nanoaod::FlatTable::IntColumn);
    
    
    globalTable->addColumn<int>("jetIdx", global_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    globalFillerList.fill(globalTable);
    
    csvTable->addColumn<int>("jetIdx", csv_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    csvFillerList.fill(csvTable);

    cpfTable->addColumn<int>("jetIdx", cpf_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    cpfFillerList.fill(cpfTable);

    npfTable->addColumn<int>("jetIdx", npf_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    npfFillerList.fill(npfTable);

    svTable->addColumn<int>("jetIdx", sv_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    svFillerList.fill(svTable);
    
    svTable_adapted->addColumn<int>("jetIdx", sv_adapted_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    svAdaptedFillerList.fill(svTable_adapted);

    muonTable->addColumn<int>("jetIdx", mu_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    muonFillerList.fill(muonTable);

    electronTable->addColumn<int>("jetIdx", elec_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    electronFillerList.fill(electronTable);
    
    iEvent.put(std::move(globalTable),     "global");
    iEvent.put(std::move(csvTable),        "csv");
    iEvent.put(std::move(cpfTable),        "cpf");
    iEvent.put(std::move(npfTable),        "npf");
    iEvent.put(std::move(svTable),         "sv");
    iEvent.put(std::move(svTable_adapted), "svAdapted");
    iEvent.put(std::move(lengthTable),     "length");
    iEvent.put(std::move(muonTable),       "muon");
    iEvent.put(std::move(electronTable),   "electron");
}

void
NANOProducer::beginStream(edm::StreamID) {
}

void
NANOProducer::endStream() {
}


void
NANOProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NANOProducer);
