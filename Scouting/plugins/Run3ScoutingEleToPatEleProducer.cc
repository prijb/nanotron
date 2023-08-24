// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingEleToPatEleProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingEleToPatEleProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingEleToPatEleProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingElectron>> electronToken_;
};

Run3ScoutingEleToPatEleProducer::Run3ScoutingEleToPatEleProducer(const edm::ParameterSet &iConfig) {
    electronToken_ = consumes<edm::View<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("eleSource"));
    produces<std::vector<pat::Electron>>();
}

Run3ScoutingEleToPatEleProducer::~Run3ScoutingEleToPatEleProducer() {}

void Run3ScoutingEleToPatEleProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingElectron>> electrons;
    iEvent.getByToken(electronToken_, electrons);

    auto patElectrons = std::make_unique<std::vector<pat::Electron>>();
    for (const auto &e: *electrons) {
        pat::Electron patE;
        math::PtEtaPhiMLorentzVectorD lv(e.pt(), e.eta(), e.phi(), e.m());
        patE.setP4(lv);
        patElectrons->push_back(patE);
    }
    iEvent.put(std::move(patElectrons));
}

void Run3ScoutingEleToPatEleProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingElectron to Pat Jet producer module");

  // input source
  iDesc.add<edm::InputTag>("eleSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingEleToPatEleProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingEleToPatEleProducer);

