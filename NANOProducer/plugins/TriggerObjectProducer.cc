// -*- C++ -*-
//
// Package:    nanotron/TriggerObjectProducer
// Class:      TriggerObjectProducer
// 
/**\class TriggerObjectProducer TriggerObjectProducer.cc nanotron/TriggerObjectProducer/plugins/TriggerObjectProducer.cc
/*
Simplified version of the MuonBPark producer with only trigger muons stored
*/
//
// Original Author:  Prijith Pradeep <prijith.babu-pradeep18@imperial.ac.uk>
//         Created:  Wed, 30 Nov 2023 10:23:30 GMT
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = false; //false;

class TriggerObjectProducer : public edm::stream::EDProducer<> {
    
public:
    
    explicit TriggerObjectProducer(const edm::ParameterSet &iConfig);
    ~TriggerObjectProducer() override {};
    virtual void produce(edm::Event&, const edm::EventSetup&);
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    
    
private:
    // const edm::ESHandle<MagneticField> bFieldHandle_; // Does not work in CMSSW_10_x (Mikael)
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;

    std::string name_;
    const double ptMin_;
    const int objId_;              //id of trigger object (+83 for muon)
    std::vector<std::string> HLTPaths_;
    
//    std::vector<std::string> L1Seeds_;
};

//Initialise producer
TriggerObjectProducer::TriggerObjectProducer(const edm::ParameterSet &iConfig):
  //bFieldHandle_(esConsumes<magField, IdealMagneticFieldRecord>()),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  magFieldToken_(esConsumes()),
  name_(iConfig.getParameter<std::string>("name")),
  ptMin_(iConfig.getParameter<double>("ptMin")),
  objId_(iConfig.getParameter<int>("objId")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths"))//,   //////////Comma
//  L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1seeds"))
{
    //Produces a table of trigger muons
    //produces<pat::MuonCollection>("triggerMuons"); 
    produces<nanoaod::FlatTable>();
    
}

//Produce trigger objects
void TriggerObjectProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    const MagneticField *magneticField_ = &iSetup.getData(magFieldToken_);

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);

    auto triggeringObjects = std::make_unique<std::vector<pat::Muon>>();
    //std::vector<pat::TriggerObjectStandAlone> triggeringObjects;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);

    std::vector<float> obj_pt, obj_eta, obj_phi, obj_mass, obj_vx, obj_vy, obj_vz;
    std::vector<int> obj_charge, obj_pdgid;
    std::vector<int> isTriggering, isTriggering_HLT_Mu4_3_LowMass, isTriggering_HLT_Mu4_LowMass_Displaced;

    //std::cout << std::endl;
    /*
    for (const std::string path: HLTPaths_){
        //This gives a generic title for the path
        char cstr[ (path+"*").size() + 1];
        strcpy( cstr, (path+"*").c_str() );    
        std::cout << "HLT Path considered: " << path.c_str() << std::endl;
    }
    */

    //Iterate over trigger objects, then iterate over paths to see if it triggers any
    int ntrig = 0;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits); 
    for(pat::TriggerObjectStandAlone obj : *triggerObjects){
        int temptrig = 0;
        int temp_hltmu4_3 = 0;
        int temp_hltmu4_disp = 0;
        //Unpack first
        obj.unpackPathNames(names);
        if(!obj.id(83)) continue; //skips non muons
        //Iterate over paths
        for (const std::string path: HLTPaths_){
            //This gives a generic title for the path
            char cstr[ (path+"*").size() + 1];
            strcpy( cstr, (path+"*").c_str() ); 
            //If the object a path, include it
            if(obj.hasPathName(cstr, true, true)){
                //std::cout << "Object passes trigger: " << cstr << std::endl;
                temptrig += 1;
                //Pick out HLT_DoubleMu4_3
                if(strcmp(cstr, "HLT_DoubleMu4_3_LowMass*")==0) temp_hltmu4_3 = 1;
                if(strcmp(cstr, "HLT_DoubleMu4_LowMass_Displaced*")==0) temp_hltmu4_disp = 1;
            }
        }
        //Now store information for objects passing any path
        if(temptrig>0){
            //std::cout << "Storing object: " <<  "  - pt: " << obj.pt() << ", eta: " << obj.eta() << ", phi: " << obj.phi() << ", mass: " << obj.mass() <<  std::endl;
            //std::cout <<  "Charge : " << obj.charge() << ", Id: " << obj.pdgId() << ", threecharge: " << obj.threeCharge() << std::endl;
            //std::cout <<  "vx : " << obj.vx() << ", vy: " << obj.vy() << ", vz: " << obj.vz() << std::endl;
            //std::cout <<  "vx : " << (obj.vertex()).x() << ", vy: " << (obj.vertex()).y() << ", vz: " << (obj.vertex()).z() << std::endl;
            //std::cout <<  "vx : " << (obj.vertex()).X() << ", vy: " << (obj.vertex()).Y() << ", vz: " << (obj.vertex()).Z() << std::endl;
            obj_pt.push_back(obj.pt());
            obj_eta.push_back(obj.eta());
            obj_phi.push_back(obj.phi());
            //Placeholder muon mass
            obj_mass.push_back(0.1057);
            obj_vx.push_back(obj.vx());
            obj_vy.push_back(obj.vy());
            obj_vz.push_back(obj.vz());
            obj_charge.push_back(obj.charge());
            obj_pdgid.push_back(obj.pdgId());

            isTriggering.push_back(1);
            isTriggering_HLT_Mu4_3_LowMass.push_back(temp_hltmu4_3);
            isTriggering_HLT_Mu4_LowMass_Displaced.push_back(temp_hltmu4_disp);
            ntrig += 1;
        }
    }

    //Now make a table
    auto tab  = std::make_unique<nanoaod::FlatTable>(ntrig, name_, false, false);
    tab->addColumn<float>("pt", obj_pt, "pt", 12);
    tab->addColumn<float>("eta", obj_eta, "eta", 12);
    tab->addColumn<float>("phi", obj_phi, "phi", 12);
    tab->addColumn<float>("vx", obj_vx, "Vertex x position", 12);
    tab->addColumn<float>("vy", obj_vy, "Vertex y position", 12);
    tab->addColumn<float>("vz", obj_vz, "Vertex z position", 12);
    tab->addColumn<int>("charge", obj_charge, "charge of HLT object");
    tab->addColumn<int>("pdgId", obj_pdgid, "ID of HLT object");
    tab->addColumn<int>("isTriggering", isTriggering, "Does it trigger any of the paths");
    tab->addColumn<int>("isTriggering_HLT_DoubleMu4_3_LowMass", isTriggering_HLT_Mu4_3_LowMass, "Triggers the low mass dimuon trigger");
    tab->addColumn<int>("isTriggering_HLT_DoubleMu4_LowMass_Displaced", isTriggering_HLT_Mu4_LowMass_Displaced, "Triggers the displaced low mass dimuon trigger");
    iEvent.put(std::move(tab));

}

void TriggerObjectProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectProducer);