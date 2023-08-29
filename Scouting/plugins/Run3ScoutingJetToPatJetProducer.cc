// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingJetToPatJetProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingJetToPatJetProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingJetToPatJetProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingPFJet>> jetsToken_;
};

Run3ScoutingJetToPatJetProducer::Run3ScoutingJetToPatJetProducer(const edm::ParameterSet &iConfig) {
    jetsToken_ = consumes<edm::View<Run3ScoutingPFJet>>(iConfig.getParameter<edm::InputTag>("jetSource"));
    produces<std::vector<pat::Jet>>("jets");
    produces<std::vector<pat::Tau>>("taus");
}

Run3ScoutingJetToPatJetProducer::~Run3ScoutingJetToPatJetProducer() {}

void Run3ScoutingJetToPatJetProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingPFJet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    auto patJets = std::make_unique<std::vector<pat::Jet>>();
    auto patTaus = std::make_unique<std::vector<pat::Tau>>();
    for (const auto &jet: *jets) {
        pat::Jet patjet;
        pat::Tau pattau;

        math::PtEtaPhiMLorentzVectorD lv(jet.pt(), jet.eta(), jet.phi(), jet.m());
        patjet.setP4(lv);
        pattau.setP4(lv);

        patJets->push_back(patjet);
        patTaus->push_back(pattau);
    }
    iEvent.put(std::move(patJets), "jets");
    iEvent.put(std::move(patTaus), "taus");
}

void Run3ScoutingJetToPatJetProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingJet to Pat Jet producer module");

  // input source
  iDesc.add<edm::InputTag>("jetSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingJetToPatJetProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingJetToPatJetProducer);

