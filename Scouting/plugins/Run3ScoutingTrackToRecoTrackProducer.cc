// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingTrackToRecoTrackProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingTrackToRecoTrackProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingTrackToRecoTrackProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingTrack>> trackToken_;
};

Run3ScoutingTrackToRecoTrackProducer::Run3ScoutingTrackToRecoTrackProducer(const edm::ParameterSet &iConfig) {
    trackToken_ = consumes<edm::View<Run3ScoutingTrack>>(iConfig.getParameter<edm::InputTag>("trackSource"));
    produces<std::vector<reco::Track>>();
}

Run3ScoutingTrackToRecoTrackProducer::~Run3ScoutingTrackToRecoTrackProducer() {}

void Run3ScoutingTrackToRecoTrackProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingTrack>> tracks;
    iEvent.getByToken(trackToken_, tracks);

    auto recoTracks = std::make_unique<std::vector<reco::Track>>();
    for (const auto &track: *tracks) {
        reco::Track::Point v(track.tk_vx(), track.tk_vy(), track.tk_vz());
        reco::Track::Vector p(
            track.tk_pt() * cos(track.tk_phi()),
            track.tk_pt() * sin(track.tk_phi()),
            track.tk_pt() * sinh(track.tk_eta())
        );
        double vec[15];
        for (auto i = 0; i < 15; i++)
            vec[i] = 1.;
        reco::TrackBase::CovarianceMatrix cov(vec, vec + 15);
        cov(0, 1) = track.tk_qoverp_lambda_cov();
        cov(0, 2) = track.tk_qoverp_phi_cov();
        cov(0, 3) = track.tk_qoverp_dxy_cov();
        cov(0, 4) = track.tk_qoverp_dsz_cov();
        cov(1, 2) = track.tk_lambda_phi_cov();
        cov(1, 3) = track.tk_lambda_dxy_cov();
        cov(1, 4) = track.tk_lambda_dsz_cov();
        cov(2, 3) = track.tk_phi_dxy_cov();
        cov(2, 4) = track.tk_phi_dsz_cov();
        cov(3, 4) = track.tk_dxy_dsz_cov();
        reco::Track tk(track.tk_chi2(), track.tk_ndof(), v, p, -1, cov);
        recoTracks->push_back(tk);
    }
    iEvent.put(std::move(recoTracks));
}

void Run3ScoutingTrackToRecoTrackProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingTrack to Track producer module");

  // input source
  iDesc.add<edm::InputTag>("trackSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingTrackToRecoTrackProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingTrackToRecoTrackProducer);

