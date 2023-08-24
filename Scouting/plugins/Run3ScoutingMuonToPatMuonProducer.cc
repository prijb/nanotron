// Run3ScoutingMuon to pat::Muon producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingMuonToPatMuonProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingMuonToPatMuonProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingMuonToPatMuonProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingMuon>> muonToken_;
};

Run3ScoutingMuonToPatMuonProducer::Run3ScoutingMuonToPatMuonProducer(const edm::ParameterSet &iConfig) {
    muonToken_ = consumes<edm::View<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonSource"));
    produces<std::vector<pat::Muon>>();
}

Run3ScoutingMuonToPatMuonProducer::~Run3ScoutingMuonToPatMuonProducer() {}

void Run3ScoutingMuonToPatMuonProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingMuon>> muons;
    iEvent.getByToken(muonToken_, muons);

    auto patMuons = std::make_unique<std::vector<pat::Muon>>();
    for (const auto &muon: *muons) {
        pat::Muon patmuon;
        math::PtEtaPhiMLorentzVectorD lv(muon.pt(), muon.eta(), muon.phi(), muon.m());
        patmuon.setP4(lv);
        patmuon.setCharge(muon.charge());
        patMuons->push_back(patmuon);
    }
    iEvent.put(std::move(patMuons));
}

void Run3ScoutingMuonToPatMuonProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingMuon to pat::Muon producer module");

  // input source
  iDesc.add<edm::InputTag>("muonSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingMuonToPatMuonProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingMuonToPatMuonProducer);

