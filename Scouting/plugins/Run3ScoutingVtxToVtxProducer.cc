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
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>
#include <TLorentzVector.h>

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
        edm::EDGetTokenT<edm::View<Run3ScoutingVertex>> pvToken_;
        edm::EDGetTokenT<edm::View<Run3ScoutingVertex>> svToken_;
        edm::EDGetTokenT<edm::View<Run3ScoutingMuon>> muonToken_;
};

Run3ScoutingVtxToVtxProducer::Run3ScoutingVtxToVtxProducer(const edm::ParameterSet &iConfig) {
    pvToken_ = consumes<edm::View<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("pvSource"));
    svToken_ = consumes<edm::View<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("svSource"));
    muonToken_ = consumes<edm::View<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonSource"));
    produces<reco::VertexCollection>("pvs");
    produces<edm::ValueMap<float>>("pvs");
    produces<std::vector<reco::VertexCompositePtrCandidate>>("svs");
}

Run3ScoutingVtxToVtxProducer::~Run3ScoutingVtxToVtxProducer() {}

void Run3ScoutingVtxToVtxProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<Run3ScoutingVertex>> pvs;
    iEvent.getByToken(pvToken_, pvs);
    edm::Handle<edm::View<Run3ScoutingVertex>> svs;
    iEvent.getByToken(svToken_, svs);
    edm::Handle<edm::View<Run3ScoutingMuon>> muons;
    iEvent.getByToken(muonToken_, muons);

    auto primaryVertices = std::make_unique<reco::VertexCollection>();
    std::vector<float> values;
    for (const auto &vertex: *pvs) {
        Point p(vertex.x(), vertex.y(), vertex.z());
        Error3 e;
        e(0, 0) = vertex.xError();
        e(1, 1) = vertex.yError();
        e(2, 2) = vertex.zError();
        reco::Vertex recoVertex(p, e, vertex.chi2(), vertex.ndof(), vertex.tracksSize());
        primaryVertices->push_back(recoVertex);
        values.push_back(1.);
    }

    auto secondaryVertices = std::make_unique<std::vector<reco::VertexCompositePtrCandidate>>();
    int iv = 0;
    for (const auto &vertex: *svs) {
        Point p(vertex.x(), vertex.y(), vertex.z());
        math::XYZTLorentzVectorD lv;
        lv.SetPxPyPzE(1., 1., 1., 1.);
        Error3 e;
        e(0, 0) = vertex.xError();
        e(1, 1) = vertex.yError();
        e(2, 2) = vertex.zError();

        TLorentzVector sVtx_Sum;
        for (const auto &muo: *muons) {
            auto *muon_iter = &muo;
            const auto& Muon_vtxIndx_inter = muon_iter->vtxIndx();
            if (std::find(Muon_vtxIndx_inter.begin(), Muon_vtxIndx_inter.end(), iv) != Muon_vtxIndx_inter.end()) {
                TLorentzVector muon_tlv;
                muon_tlv.SetPtEtaPhiM(muon_iter->pt(), muon_iter->eta(), muon_iter->phi(), muon_iter->m());
                sVtx_Sum += muon_tlv;
            }
        }
        reco::VertexCompositePtrCandidate recoVertex(1, lv, p, e, vertex.chi2(), vertex.ndof());
        recoVertex.setMass(sVtx_Sum.M());
        secondaryVertices->push_back(recoVertex);
        iv++;
    }

    // insert PVs with a dummy score
    auto pvcol = iEvent.put(std::move(primaryVertices), "pvs");
    auto vertexScoreOutput = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler vertexScoreFiller(*vertexScoreOutput);
    vertexScoreFiller.insert(pvcol, values.begin(), values.end());
    vertexScoreFiller.fill();
    iEvent.put(std::move(vertexScoreOutput), "pvs");

    // insert SVs
    iEvent.put(std::move(secondaryVertices), "svs");
}

void Run3ScoutingVtxToVtxProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("Run3ScoutingJet to Pat Jet producer module");

  // input source
  iDesc.add<edm::InputTag>("pvSource", edm::InputTag("no default"))->setComment("input collection");
  iDesc.add<edm::InputTag>("svSource", edm::InputTag("no default"))->setComment("input collection");
  iDesc.add<edm::InputTag>("muonSource", edm::InputTag("no default"))->setComment("input collection");

  descriptions.add("Run3ScoutingVtxToVtxProducer", iDesc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingVtxToVtxProducer);

