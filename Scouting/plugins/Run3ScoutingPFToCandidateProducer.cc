// Run3ScoutingMuon to pat::Muon producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <math.h>
#include <unordered_map>
#include <iostream>
#include <regex>

typedef math::XYZPoint Point;

class Run3ScoutingPFToCandidateProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingPFToCandidateProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingPFToCandidateProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> particleToken_;
        edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexToken_;
};

Run3ScoutingPFToCandidateProducer::Run3ScoutingPFToCandidateProducer(const edm::ParameterSet &iConfig) {
    particleToken_ = consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("particleSource"));
    vertexToken_ = consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("vertexSource"));
    produces<std::vector<reco::CompositeCandidate>>();
}

Run3ScoutingPFToCandidateProducer::~Run3ScoutingPFToCandidateProducer() {}

void Run3ScoutingPFToCandidateProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<std::vector<Run3ScoutingParticle>> particles;
    iEvent.getByToken(particleToken_, particles);
    edm::Handle<std::vector<Run3ScoutingVertex>> vertices;
    iEvent.getByToken(vertexToken_, vertices);

    auto candidates = std::make_unique<std::vector<reco::CompositeCandidate>>();

    for (const auto &part: *particles) {
        if (part.vertex() >= 0) {
          auto vertex = vertices->at(part.vertex());
          Point p(vertex.x(), vertex.y(), vertex.z());
          math::PtEtaPhiMLorentzVectorD lv(part.pt(), part.eta(), part.phi(), 1);
          reco::CompositeCandidate cand(0, lv, p, part.pdgId());
          candidates->push_back(cand);
        } else {
          Point p(-1, -1, -1);
          math::PtEtaPhiMLorentzVectorD lv(part.pt(), part.eta(), part.phi(), 1);
          reco::CompositeCandidate cand(0, lv, p, part.pdgId());
          candidates->push_back(cand);
        }
    }
    iEvent.put(std::move(candidates));
}

void Run3ScoutingPFToCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingParticle to reco::Candidate producer module");

  // input source
  iDesc.add<edm::InputTag>("particleSource", edm::InputTag("no default"))->setComment("input collection");
  iDesc.add<edm::InputTag>("vertexSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingPFToCandidateProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingPFToCandidateProducer);

