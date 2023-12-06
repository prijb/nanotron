// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

class Run3ScoutingPhotonToPatPhotonProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingPhotonToPatPhotonProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingPhotonToPatPhotonProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingPhoton>> photonToken_;
};

Run3ScoutingPhotonToPatPhotonProducer::Run3ScoutingPhotonToPatPhotonProducer(const edm::ParameterSet &iConfig) {
    photonToken_ = consumes<edm::View<Run3ScoutingPhoton>>(iConfig.getParameter<edm::InputTag>("photonSource"));
    produces<std::vector<pat::Photon>>();
}

Run3ScoutingPhotonToPatPhotonProducer::~Run3ScoutingPhotonToPatPhotonProducer() {}

void Run3ScoutingPhotonToPatPhotonProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingPhoton>> photons;
    iEvent.getByToken(photonToken_, photons);

    auto patPhotons = std::make_unique<std::vector<pat::Photon>>();
    for (const auto &p: *photons) {
        pat::Photon patP;
        math::PtEtaPhiMLorentzVectorD lv(p.pt(), p.eta(), p.phi(), p.m());
        patP.setP4(lv);
        patPhotons->push_back(patP);
    }
    iEvent.put(std::move(patPhotons));
}

void Run3ScoutingPhotonToPatPhotonProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingPhoton to pat::Photon producer module");

  // input source
  iDesc.add<edm::InputTag>("photonSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingPhotonToPatPhotonProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingPhotonToPatPhotonProducer);

