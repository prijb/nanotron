//
//
//
//

#include "DataFormats/Common/interface/Wrapper.h"

#include "nanotron/DataFormats/interface/OnionTagFeatures.h"
#include "nanotron/DataFormats/interface/OnionTagInfo.h"
#include "nanotron/DataFormats/interface/DisplacedGenVertex.h"

#include "nanotron/DataFormats/interface/MCLabel.h"
#include "nanotron/DataFormats/interface/MCGenDecayInfo.h"
#include "nanotron/DataFormats/interface/MCGhostFlavourInfo.h"


namespace {

    struct dictionary {
        
        std::vector<reco::FeaturesTagInfo<nanotron::OnionTagFeatures>> dummy0;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<nanotron::OnionTagFeatures>>> dummy1;
        
        reco::FeaturesTagInfo<nanotron::OnionTagFeatures> dummy2;
        edm::Wrapper<reco::FeaturesTagInfo<nanotron::OnionTagFeatures>> dummy3;
        
        nanotron::OnionTagFeatures dummy4;
        nanotron::JetFeatures dummy5;
        nanotron::SecondaryVertexFeatures dummy6;
        nanotron::ChargedCandidateFeatures dummy7;
        nanotron::NeutralCandidateFeatures dummy8;
        nanotron::DisplacedGenVertexCollection dummy9;


        edm::Wrapper<nanotron::DisplacedGenVertexCollection> dummy10;
        
        edm::Ptr<nanotron::DisplacedGenVertex> dummy11;
        edm::Wrapper<edm::Ptr<nanotron::DisplacedGenVertex>> dummy12;
        
        edm::Ptr<nanotron::DisplacedGenVertexCollection> dummy13;
        edm::Wrapper<edm::Ptr<nanotron::DisplacedGenVertexCollection>> dummy14;

        edm::PtrVector<nanotron::DisplacedGenVertexCollection> dummy15;
        edm::Wrapper<edm::PtrVector<nanotron::DisplacedGenVertexCollection>> dummy16;
        
        edm::PtrVector<reco::GenParticle> dummy17;
        edm::Wrapper<edm::PtrVector<reco::GenParticle>> dummy18;
        
        /*
        nanotron::MCLabel dummy19;
        reco::FeaturesTagInfo<nanotron::MCLabel> dummy20;
        
        std::vector<reco::FeaturesTagInfo<nanotron::MCLabel>> dummy21;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<nanotron::MCLabel>>> dummy22;
        
        nanotron::MCGenDecayInfo dummy23;
        std::vector<nanotron::MCGenDecayInfo> dummy24;
        edm::Wrapper<std::vector<nanotron::MCGenDecayInfo>> dummy25;
        
        nanotron::MCGhostFlavourInfo dummy26;
        std::vector<nanotron::MCGhostFlavourInfo> dummy27;
        edm::Wrapper<std::vector<nanotron::MCGhostFlavourInfo>> dummy28;
        edm::ValueMap<nanotron::MCGhostFlavourInfo> dummy29;
        edm::Wrapper<edm::ValueMap<nanotron::MCGhostFlavourInfo>> dummy30;
        */

        nanotron::ElectronCandidateFeatures dummy19;
        nanotron::MuonCandidateFeatures dummy20;
    };
}
