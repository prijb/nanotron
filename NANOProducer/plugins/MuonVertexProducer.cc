// -*- C++ -*-
//
// Package:    nanotron/MuonVertexProducer
// Class:      MuonVertexProducer
// 
/**\class MuonVertexProducer MuonVertexProducer.cc nanotron/MuonVertexProducer/plugins/MuonVertexProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Tue, 24 Jan 2023 23:20:36 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TMath.h"
#include <Math/Vector4D.h>


#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

//
// class declaration
//

class MuonVertexProducer : public edm::stream::EDProducer<> {
   public:
      explicit MuonVertexProducer(const edm::ParameterSet&);
      ~MuonVertexProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;
      
      edm::InputTag _muonInputTag;
      edm::EDGetTokenT<std::vector<pat::Muon>> _muonToken;
      const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbESToken_;
      const std::string  svName_;
      const double ptMin_,dlenMin_,dlenSigMin_;
      const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
      const StringCutObjectSelector<reco::Vertex> svCut_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
MuonVertexProducer::MuonVertexProducer(const edm::ParameterSet& iConfig):
    _muonInputTag(iConfig.getParameter<edm::InputTag>("srcMuon")),
    _muonToken(consumes<std::vector<pat::Muon>>(_muonInputTag)),
    ttbESToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    svName_(iConfig.getParameter<std::string>("svName") ),
    ptMin_(iConfig.getParameter<double>("ptMin") ),
    dlenMin_(iConfig.getParameter<double>("dlenMin") ),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin") ),
    pvs_(consumes<std::vector<reco::Vertex>>( iConfig.getParameter<edm::InputTag>("pvSrc") )),
    svCut_(iConfig.getParameter<std::string>("svCut") , true)
{
   produces<std::vector<reco::Vertex>>();
   produces<nanoaod::FlatTable>("svs");
}


MuonVertexProducer::~MuonVertexProducer()
{

   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
       // Pick pair of muons with smallest vertex chi square fit for all collection combos
    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(_muonToken, muons);
    std::vector<reco::TrackRef> muTracks{};
    std::vector<pat::Muon> muObjs{};

    edm::Handle<std::vector<reco::Vertex>> pvsIn;
    iEvent.getByToken(pvs_, pvsIn);

    std::vector<float> dlen,dlenSig,pAngle,dxy,dxySig,x,y,z,ndof,chi2,origMass,propMass,mu1pt,mu2pt,mu1phi,mu2phi,mu1eta,mu2eta;
    std::vector<int> mu1index,mu2index,charge;
    VertexDistance3D vdist;
    VertexDistanceXY vdistXY;

    size_t nGoodSV=0;
    float muonMass=0.1057;
    const auto & PV0 = pvsIn->front();
    std::vector<int> origIndex;
    for (size_t i = 0; i < muons->size(); i++) {
        reco::TrackRef track_i = muons->at(i).muonBestTrack();
        if (track_i.isNonnull() && track_i->pt() > ptMin_) {
            muTracks.emplace_back(track_i);
            muObjs.emplace_back(muons->at(i));
            origIndex.emplace_back(i);
        }
    }

    std::sort(origIndex.begin(), origIndex.end(),
            [&muTracks](const int& a, const int& b)
            {
            return muTracks[a]->pt() > muTracks[b]->pt();
            });
    sort(muTracks.begin(), muTracks.end(), [](const auto & l, const auto & r) {
            return l->pt() > r->pt();
            });
    sort(muObjs.begin(), muObjs.end(), [](const auto & l, const auto & r) {
            reco::TrackRef lt = l.muonBestTrack();
            reco::TrackRef rt = r.muonBestTrack();
            return lt->pt() > rt->pt();
            });

    // edm::ESHandle<TransientTrackBuilder> theB;
    // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
    const TransientTrackBuilder* theB = &iSetup.getData(ttbESToken_);
    KalmanVertexFitter kvf(true);
    auto vertices = std::make_unique<std::vector<reco::Vertex>>();

    for (size_t i = 0; i < muTracks.size(); i++) {
        for (size_t j = i+1; j < muTracks.size(); j++) {
            reco::TrackRef muon_i, muon_j;
            if (i < muTracks.size())
                muon_i = muTracks[i];
            if (j < muTracks.size())
                muon_j = muTracks[j];

            TransientVertex tv;
            if (muon_i.isNonnull() && muon_j.isNonnull() && i != j) {
            //Removed OS requirement
            //if (muon_i.isNonnull() && muon_j.isNonnull() && i != j && muObjs[i].charge() != muObjs[j].charge()) {
                std::vector<reco::TransientTrack> transient_tracks{};
                transient_tracks.push_back(theB->build(muon_i));
                transient_tracks.push_back(theB->build(muon_j));
                tv = kvf.vertex(transient_tracks);

                // float vxy = -9999;
                // float sigma_vxy = -9999;
                // float vtx_chi2 = 999999;
                // float vz = -9999;
                // float dr = -9999;

                if (tv.isValid()){
                    reco::Vertex vertex = reco::Vertex(tv);
                    if (svCut_(vertex)) {
                        Measurement1D dl= vdist.distance(PV0,VertexState(RecoVertex::convertPos(vertex.position()),RecoVertex::convertError(vertex.error())));
                        if(dl.value() > dlenMin_ and dl.significance() > dlenSigMin_){
                            double dx = (-PV0.x() + vertex.x()), dy = (-PV0.y() + vertex.y()), dz = (-PV0.z() + vertex.z());
                            double pmag = sqrt(vertex.p4(muonMass).px()*vertex.p4(muonMass).px()+vertex.p4(muonMass).py()*vertex.p4(muonMass).py()+vertex.p4(muonMass).pz()*vertex.p4(muonMass).pz());
                            double pdotv = (dx * vertex.p4(muonMass).px() + dy*vertex.p4(muonMass).py() + dz*vertex.p4(muonMass).pz())/(pmag*sqrt(dx*dx + dy*dy + dz*dz));
                            // std::cout << vertex.p4(0.107).mass() << std::endl;
                            // vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
                            // sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                            //         vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
                            // //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
                            // vtx_chi2 = vertex.normalizedChi2();
                            // vz = vertex.z();
                            // dr = reco::deltaR(*muon_i, *muon_j);
                            ROOT::Math::PtEtaPhiMVector propagatedP4(0,0,0,0);
                            for(auto trans : transient_tracks){
                                GlobalPoint vert(vertex.x(),vertex.y(),vertex.z());
                                TrajectoryStateClosestToPoint  traj = trans.trajectoryStateClosestToPoint(vert);
                                GlobalVector momentumPropagated = traj.momentum();
                                // std::cout << momentumPropagated << std::endl;
                                ROOT::Math::PtEtaPhiMVector muTrans = ROOT::Math::PtEtaPhiMVector(momentumPropagated.perp(),momentumPropagated.eta(),momentumPropagated.phi(),muonMass);
                                propagatedP4 += muTrans;
                            } 
                            vertices->push_back(vertex);
                            dlen.push_back(dl.value());
                            dlenSig.push_back(dl.significance());
                            origMass.push_back(vertex.p4(muonMass).mass());
                            propMass.push_back(propagatedP4.mag());
                            pAngle.push_back(std::acos(pdotv));
                            Measurement1D d2d = vdistXY.distance(PV0, VertexState(RecoVertex::convertPos(vertex.position()), RecoVertex::convertError(vertex.error())));
                            dxy.push_back(d2d.value());
                            dxySig.push_back(d2d.significance());
                            x.push_back(vertex.x());
                            y.push_back(vertex.y());
                            z.push_back(vertex.z());
                            ndof.push_back(vertex.ndof());
                            chi2.push_back(vertex.normalizedChi2());
                            mu1pt.push_back(muObjs[i].pt());
                            mu2pt.push_back(muObjs[j].pt());
                            mu1phi.push_back(muObjs[i].phi());
                            mu2phi.push_back(muObjs[j].phi());
                            mu1eta.push_back(muObjs[i].eta());
                            mu2eta.push_back(muObjs[j].eta());
                            mu1index.push_back(origIndex[i]);
                            mu2index.push_back(origIndex[j]);
                            charge.push_back(muObjs[i].charge() + muObjs[j].charge());
                            // std::cout << muObjs[i]->pt() << " " << muons->at(mu1index[-1])->pt() << std::endl;
                            nGoodSV++;
                        }
                    }
                }
            }
        }
    }
    //
    // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer 
    auto svsTable = std::make_unique<nanoaod::FlatTable>(nGoodSV,svName_,false);
    svsTable->addColumn<float>("dlen",dlen,"decay length in cm",23);
    svsTable->addColumn<float>("dlenSig",dlenSig,"decay length significance", 20);
    svsTable->addColumn<float>("dxy", dxy, "2D decay length in cm", 20);
    svsTable->addColumn<float>("dxySig", dxySig, "2D decay length significance", 20);
    svsTable->addColumn<float>("x",x,  "secondary vertex X position, in cm",23);
    svsTable->addColumn<float>("y",y,  "secondary vertex Y position, in cm",23);
    svsTable->addColumn<float>("z",z,  "secondary vertex Z position, in cm",23);
    svsTable->addColumn<float>("ndof",ndof,"number of degrees of freedom",10);
    svsTable->addColumn<float>("chi2",chi2, "reduced chi2, i.e. chi/ndof",10);
    svsTable->addColumn<float>("pAngle",pAngle,"pointing angle, i.e. acos(p_SV * (SV - PV)) ",23);
    svsTable->addColumn<float>("origMass",origMass,"original mass from the vertex p4",23);
    svsTable->addColumn<float>("mass",propMass,"mass propagated to the vertex position",23);
    svsTable->addColumn<float>("mu1pt",mu1pt,  "lead muon pt for vertex",23);
    svsTable->addColumn<float>("mu2pt",mu2pt,  "second muon pt for vertex",23);
    svsTable->addColumn<float>("mu1phi",mu1phi,  "lead muon phi for vertex",23);
    svsTable->addColumn<float>("mu2phi",mu2phi,  "second muon phi for vertex",23);
    svsTable->addColumn<float>("mu1eta",mu1eta,  "lead muon eta for vertex",23);
    svsTable->addColumn<float>("mu2eta",mu2eta,  "second muon eta for vertex",23);
    svsTable->addColumn<int>("mu1index",mu1index,  "lead muon index for vertex");
    svsTable->addColumn<int>("mu2index",mu2index,  "second muon index for vertex");
    svsTable->addColumn<int>("charge",charge,  "vertex charge");


    iEvent.put(std::move(vertices));
    iEvent.put(std::move(svsTable),"svs");
    /* This is an event example
    //Read 'ExampleData' from the Event
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);

    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    iEvent.put(std::make_unique<ExampleData2>(*pIn));
    */

    /* this is an EventSetup example
    //Read SetupData from the SetupRecord in the EventSetup
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    */

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
    void
MuonVertexProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuonVertexProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   MuonVertexProducer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   MuonVertexProducer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   MuonVertexProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   MuonVertexProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonVertexProducer);
