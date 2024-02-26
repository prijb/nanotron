// -*- C++ -*-
//
// Package:    nanotron/Run3ScoutingVtxTableProducer
// Class:      Run3ScoutingVtxTableProducer
// 
/**\class Run3ScoutingVtxTableProducer Run3ScoutingVtxTableProducer.cc nanotron/Run3ScoutingVtxTableProducer/plugins/Run3ScoutingVtxTableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/


// system include files
#include <memory>

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

//Flat table and vertex distance calculation utilities
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>
#include <TLorentzVector.h>

typedef math::XYZPoint Point;
typedef math::Error<3>::type Error3;

//
// class declaration
//

//Note: I did it like this because reco::VertexCompositePtrCandidate doesn't like adding user ints and floats?
//Unlike the simple candidate flat table producer used in the CMSSW repo, this also includes muon cross reference values, vertex
//dxy, dxySig, dlen, dlenSig and pValue

class Run3ScoutingVtxTableProducer : public edm::stream::EDProducer<> {
   public:
        Run3ScoutingVtxTableProducer(const edm::ParameterSet &iConfig);
        ~Run3ScoutingVtxTableProducer() override;
        void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

        static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    private:
        edm::EDGetTokenT<edm::View<Run3ScoutingVertex>> pvToken_;
        edm::EDGetTokenT<edm::View<Run3ScoutingVertex>> svToken_;
        edm::EDGetTokenT<edm::View<Run3ScoutingMuon>> muonToken_;
        const std::string svName_;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
Run3ScoutingVtxTableProducer::Run3ScoutingVtxTableProducer(const edm::ParameterSet& iConfig):
    pvToken_(consumes<edm::View<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("pvSource"))),
    svToken_(consumes<edm::View<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("svSource"))),
    muonToken_(consumes<edm::View<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonSource"))),
    svName_(iConfig.getParameter<std::string>("svName") )
{
    produces<nanoaod::FlatTable>("svs");
}

Run3ScoutingVtxTableProducer::~Run3ScoutingVtxTableProducer(){}

void
Run3ScoutingVtxTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<Run3ScoutingVertex>> pvs;
    iEvent.getByToken(pvToken_, pvs);
    edm::Handle<edm::View<Run3ScoutingVertex>> svs;
    iEvent.getByToken(svToken_, svs);
    edm::Handle<edm::View<Run3ScoutingMuon>> muons;
    iEvent.getByToken(muonToken_, muons);

    //SV columns to be populated
    std::vector<float> dlen,dlenSig,pAngle,dxy,dxySig,x,y,z,xError,yError,zError,chi2,ndof,mass;
    std::vector<int> trksize,mu1idx,mu2idx,mu3idx,mu4idx,nMuonMatched,isValidVtx;
    //Vertex distance calculation
    VertexDistance3D vdist;
    VertexDistanceXY vdistXY;

    //Produce primary vertices (first one will be used for SV dxy and dxySig)
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

    //Get the first PV
    const auto& PV0 = primaryVertices->front();
    
    //Produce a secondary vertex composite candidate to get dlen, dlenSig, pAngle, dxy and dxySig
    int iv=0;
    for(const auto &vertex: *svs){
        auto secondaryVertex = std::make_unique<reco::VertexCompositePtrCandidate>();
        Point p(vertex.x(), vertex.y(), vertex.z());
        math::XYZTLorentzVectorD lv;
        lv.SetPxPyPzE(0., 0., 0., 0.);
        Error3 e;
        e(0, 0) = vertex.xError();
        e(1, 1) = vertex.yError();
        e(2, 2) = vertex.zError();

        //Vector that gets the index reference to at most, four muons
        std::vector<int> muidx_vec(4,-1);
        int nMuonMatch = 0;
        int muidx = 0;
        //Placeholder value for mass
        float vertex_mass = -1;

        //Iterate over muons in event to get cross reference and mass of SV
        TLorentzVector sVtx_Sum;
        for (const auto &muo: *muons) {
            auto *muon_iter = &muo;
            const auto& Muon_vtxIndx_inter = muon_iter->vtxIndx();

            if (std::find(Muon_vtxIndx_inter.begin(), Muon_vtxIndx_inter.end(), iv) != Muon_vtxIndx_inter.end()) {
                //Muon is found
                if(nMuonMatch<4) muidx_vec[nMuonMatch] = muidx;
                TLorentzVector muon_tlv;
                //Putting this because idk if it's safe to use TLorentzVector
                math::XYZTLorentzVectorD muon_tlvd; 
                muon_tlv.SetPtEtaPhiM(muon_iter->pt(), muon_iter->eta(), muon_iter->phi(), muon_iter->m());
                muon_tlvd.SetPxPyPzE(muon_tlv.Px(), muon_tlv.Py(), muon_tlv.Pz(), muon_tlv.E());
                sVtx_Sum += muon_tlv;
                lv += muon_tlvd;
                nMuonMatch += 1;
            }
            muidx += 1;
        }
        if(nMuonMatch>0) vertex_mass = sVtx_Sum.M();
        reco::VertexCompositePtrCandidate recoVertex(1, lv, p, e, vertex.chi2(), vertex.ndof());

        //Now measure distances from PV (recoVertex is just vertex as a CompositePtrCandidate)
        Measurement1D dl= vdist.distance(PV0,VertexState(RecoVertex::convertPos(recoVertex.position()),RecoVertex::convertError(recoVertex.error())));
        Measurement1D d2d = vdistXY.distance(PV0, VertexState(RecoVertex::convertPos(recoVertex.position()), RecoVertex::convertError(recoVertex.error())));
        double dx = (-PV0.x() + vertex.x()), dy = (-PV0.y() + vertex.y()), dz = (-PV0.z() + vertex.z());
        double pmag = sqrt(lv.px()*lv.px()+lv.py()*lv.py()+lv.pz()*lv.pz());
        double pdotv = (dx * lv.px() + dy*lv.py() + dz*lv.pz())/(pmag*sqrt(dx*dx + dy*dy + dz*dz));
        
        //SV quantities
        dlen.push_back(dl.value());
        dlenSig.push_back(dl.significance());
        pAngle.push_back(std::acos(pdotv));
        x.push_back(vertex.x());
        y.push_back(vertex.y());
        z.push_back(vertex.z());
        dxy.push_back(d2d.value());
        dxySig.push_back(d2d.significance());
        xError.push_back(vertex.xError());
        yError.push_back(vertex.yError());
        zError.push_back(vertex.zError());
        chi2.push_back(vertex.chi2());
        mass.push_back(vertex_mass);
        ndof.push_back(vertex.ndof());
        //Ints
        trksize.push_back(vertex.tracksSize());
        mu1idx.push_back(muidx_vec[0]);
        mu2idx.push_back(muidx_vec[1]);
        mu3idx.push_back(muidx_vec[2]);
        mu4idx.push_back(muidx_vec[3]);
        nMuonMatched.push_back(nMuonMatch);
        isValidVtx.push_back(int(vertex.isValidVtx()));

        iv++;
    }
    

    // Make the SV table
    auto svsTable = std::make_unique<nanoaod::FlatTable>(iv,svName_,false);
    svsTable->addColumn<float>("dlen",dlen,"decay length in cm",23);
    svsTable->addColumn<float>("dlenSig",dlenSig,"decay length significance", 20);
    svsTable->addColumn<float>("dxy", dxy, "2D decay length in cm", 20);
    svsTable->addColumn<float>("dxySig", dxySig, "2D decay length significance", 20);
    svsTable->addColumn<float>("x",x,  "secondary vertex X position, in cm",23);
    svsTable->addColumn<float>("y",y,  "secondary vertex Y position, in cm",23);
    svsTable->addColumn<float>("z",z,  "secondary vertex Z position, in cm",23);
    svsTable->addColumn<float>("xError",xError,  "secondary vertex X position error, in cm",23);
    svsTable->addColumn<float>("yError",yError,  "secondary vertex Y position error, in cm",23);
    svsTable->addColumn<float>("zError",zError,  "secondary vertex Z position error, in cm",23);
    svsTable->addColumn<float>("ndof",ndof,"number of degrees of freedom",10);
    svsTable->addColumn<float>("chi2",chi2, "reduced chi2, i.e. chi/ndof",10);
    svsTable->addColumn<float>("pAngle",pAngle,"pointing angle, i.e. acos(p_SV * (SV - PV)) ",23);
    svsTable->addColumn<float>("mass",mass,"mass from vertex matched muons",23);
    svsTable->addColumn<int>("mu1idx",mu1idx,  "lead muon index for vertex");
    svsTable->addColumn<int>("mu2idx",mu2idx,  "second muon index for vertex");
    svsTable->addColumn<int>("mu3idx",mu3idx,  "third muon index for vertex");
    svsTable->addColumn<int>("mu4idx",mu4idx,  "fourth muon index for vertex");
    svsTable->addColumn<int>("nMuonMatched",nMuonMatched,  "Number of matched muons");
    svsTable->addColumn<int>("isValidVtx",isValidVtx,  "Is a valid vertex");

    iEvent.put(std::move(svsTable),"svs");

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Run3ScoutingVtxTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Run3ScoutingVtxTableProducer);
