// Run3ScoutingPFJet to pat::Jet producer

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

typedef math::XYZPoint Point;
typedef math::Error<3>::type Error3;

class Run3ScoutingVtxToVtxProducer: public edm::stream::EDProducer<>
{
    public:
        Run3ScoutingVtxToVtxProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingVtxToVtxProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingVertex>> vertexToken_;
};

Run3ScoutingVtxToVtxProducer::Run3ScoutingVtxToVtxProducer(const edm::ParameterSet &iConfig) {
    vertexToken_ = consumes<edm::View<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("vertexSource"));
    produces<reco::VertexCollection>();
}

Run3ScoutingVtxToVtxProducer::~Run3ScoutingVtxToVtxProducer() {}

void Run3ScoutingVtxToVtxProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingVertex>> vertices;
    iEvent.getByToken(vertexToken_, vertices);

    auto recoVertices = std::make_unique<reco::VertexCollection>();
    for (const auto &vertex: *vertices) {
        Point p(vertex.x(), vertex.y(), vertex.z());
        Error3 e;
        reco::Vertex recoVertex(p, e, vertex.chi2(), vertex.ndof(), vertex.tracksSize());
        recoVertices->push_back(recoVertex);
    }
    iEvent.put(std::move(recoVertices));
}

void Run3ScoutingVtxToVtxProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingJet to Pat Jet producer module");

  // input source
  iDesc.add<edm::InputTag>("vertexSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingVtxToVtxProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingVtxToVtxProducer);

