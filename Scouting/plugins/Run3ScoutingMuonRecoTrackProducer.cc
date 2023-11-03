// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingMuonRecoTrackProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingMuonRecoTrackProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingMuonRecoTrackProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonToken_;
};

Run3ScoutingMuonRecoTrackProducer::Run3ScoutingMuonRecoTrackProducer(const edm::ParameterSet &iConfig) {
    muonToken_ = consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonSource"));
    produces<std::vector<reco::Track>>();
}

Run3ScoutingMuonRecoTrackProducer::~Run3ScoutingMuonRecoTrackProducer() {}

void Run3ScoutingMuonRecoTrackProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<std::vector<Run3ScoutingMuon>> muons;
    iEvent.getByToken(muonToken_, muons);

    auto recoTracks = std::make_unique<std::vector<reco::Track>>();
    for (const auto &muon: *muons) {
        reco::Track::Point v(muon.trk_vx(), muon.trk_vy(), muon.trk_vz());
        reco::Track::Vector p(
            muon.trk_pt() * cos(muon.trk_phi()),
            muon.trk_pt() * sin(muon.trk_phi()),
            muon.trk_pt() * sinh(muon.trk_eta())
        );
        double vec[15];
        for (auto i = 0; i < 15; i++)
            vec[i] = 1.;
        reco::TrackBase::CovarianceMatrix cov(vec, vec + 15);
        cov(0, 0) = std::pow(muon.trk_qoverpError(),2);
        cov(0, 1) = muon.trk_qoverp_lambda_cov();
        cov(0, 2) = muon.trk_qoverp_phi_cov();
        cov(0, 3) = muon.trk_qoverp_dxy_cov();
        cov(0, 4) = muon.trk_qoverp_dsz_cov();
        cov(1, 1) = std::pow(muon.trk_lambdaError(),2);
        cov(1, 2) = muon.trk_lambda_phi_cov();
        cov(1, 3) = muon.trk_lambda_dxy_cov();
        cov(1, 4) = muon.trk_lambda_dsz_cov();
        cov(2, 2) = std::pow(muon.trk_phiError(),2);
        cov(2, 3) = muon.trk_phi_dxy_cov();
        cov(2, 4) = muon.trk_phi_dsz_cov();
        cov(3, 3) = std::pow(muon.trk_dxyError(),2);
        cov(3, 4) = muon.trk_dxy_dsz_cov();
        cov(4, 4) = std::pow(muon.trk_dszError(),2);
        reco::Track trk(muon.trk_chi2(), muon.trk_ndof(), v, p, muon.charge(), cov);

        for (auto &hit : muon.trk_hitPattern().hitPattern) {
            switch (hit % 4) {
                case 0:
                    trk.appendHitPattern(hit, TrackingRecHit::valid);
                    break;
                case 1:
                    trk.appendHitPattern(hit, TrackingRecHit::missing);
                    break;
                case 2:
                    trk.appendHitPattern(hit, TrackingRecHit::inactive);
                    break;
                case 3:
                    trk.appendHitPattern(hit, TrackingRecHit::bad);
                    break;
            }
        }
        recoTracks->push_back(trk);
    }
    iEvent.put(std::move(recoTracks));
}

void Run3ScoutingMuonRecoTrackProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingTrack to Track producer module");

  // input source
  iDesc.add<edm::InputTag>("muonSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingMuonRecoTrackProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingMuonRecoTrackProducer);

