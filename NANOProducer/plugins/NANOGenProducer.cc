//
//
//
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "nanotron/DataFormats/interface/OnionTagInfo.h"
#include "nanotron/DataFormats/interface/MCLabel.h"
#include "nanotron/DataFormats/interface/MCLabelInfo.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <Math/Vector4D.h>

#include "FlatTableFiller.h"

//
// class declaration
//

#define LABEL(name) \
    std::pair<std::string, nanotron::MCLabel::Type>( #name , nanotron::MCLabel::Type:: name ) 

/*  
#define TAUDECAY(name) \
    std::pair<std::string, nanotron::MCLabel::TauDecay>( #name , nanotron::MCLabel::TauDecay:: name ) 
*/

class NANOGenProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOGenProducer(const edm::ParameterSet&);
      ~NANOGenProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      std::vector<std::pair<std::string,nanotron::MCLabel::Type>> labelTypes;
      //std::vector<std::pair<std::string,nanotron::MCLabel::TauDecay>> tauDecayTypes;
      
      PropertyList<nanotron::MCLabel> labelProperties;

   private:
      const edm::EDGetTokenT<edm::View<pat::Jet>> _jet_src;
      const edm::EDGetTokenT<std::vector<reco::MCLabelInfo>> _label_src;

      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOGenProducer::NANOGenProducer(const edm::ParameterSet& iConfig):
    _jet_src(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))), 
    _label_src(consumes<std::vector<reco::MCLabelInfo>>(iConfig.getParameter<edm::InputTag>("srcLabels")))
{

    labelTypes = {
        LABEL(isPU),
        LABEL(isB),
        LABEL(isBB),
        LABEL(isLeptonic_B),
        LABEL(isC),
        LABEL(isCC),
        LABEL(isLeptonic_C),
        LABEL(isS),
        LABEL(isUD),
        LABEL(isG),
        LABEL(isPrompt_MU),
        LABEL(isPrompt_E),
        LABEL(isPrompt_PHOTON),
        LABEL(isPrompt_TAU),

        LABEL(isMC_RAD), //no flavour match (likely from wide angle radiation)
        LABEL(isMC_MU), //prompt lepton
        LABEL(isMC_E), //prompt lepton
        LABEL(isMC_Q), //single light quark
        LABEL(isMC_QMU), //single light quark + prompt lepton
        LABEL(isMC_QE), //single light quark + prompt lepton
        LABEL(isMC_QQ), //double light quark
        LABEL(isMC_QQMU), //double light quark + prompt lepton
        LABEL(isMC_QQE), //double light quark + prompt lepton
        LABEL(isMC_B), //single b/c quark
        LABEL(isMC_BMU), //single b/c quark + prompt lepton
        LABEL(isMC_BE), //single b/c quark + prompt lepton
        LABEL(isMC_BB), //double b/c quark
        LABEL(isMC_BBMU), //double b/c quark + prompt lepton
        LABEL(isMC_BBE), //double b/c quark + prompt lepton
        LABEL(isMC_TAU),
        LABEL(isMC_QTAU),
        LABEL(isMC_QQTAU),
        LABEL(isMC_BTAU),
        LABEL(isMC_BBTAU),

        LABEL(isMC_PHOTON),
        LABEL(isMC_QPHOTON),
        LABEL(isMC_QQPHOTON),
        LABEL(isMC_BPHOTON),
        LABEL(isMC_BBPHOTON),
        LABEL(isUndefined)
    };
    
    
    /*
    tauDecayTypes = {
        TAUDECAY(NO_TAU),     // no tau decay
        TAUDECAY(INVISIBLE),  // tau decay but not reconstructable
        TAUDECAY(E),          // to electron
        TAUDECAY(MU),         // to muon
        TAUDECAY(H),          // 1 charged hadron
        TAUDECAY(H_1PI0),     // 1 charged hadron + pi0(->2gamma)
        TAUDECAY(H_XPI0),     // 1 charged hadron + 2 or more pi0(->2gamma)
        TAUDECAY(HHH),        // 3 charged hadrons
        TAUDECAY(HHH_XPI0)    // 3 charged hadron + 1 or more pi0(->2gamma)
    };
    */

    
    labelProperties = {
        PROPERTY(nanotron::MCLabel, partonFlavor, "doc"),
        PROPERTY(nanotron::MCLabel, hadronFlavor, "doc"),
        PROPERTY(nanotron::MCLabel, llpId,        "doc"),
        PROPERTY(nanotron::MCLabel, llp_mass,     "doc"),
        PROPERTY(nanotron::MCLabel, llp_pt,       "doc"),

        PROPERTY(nanotron::MCLabel, displacement,    "doc"),
        PROPERTY(nanotron::MCLabel, displacement_xy, "doc"),
        PROPERTY(nanotron::MCLabel, displacement_z,  "doc"),
        PROPERTY(nanotron::MCLabel, decay_angle,     "doc"),
        PROPERTY(nanotron::MCLabel, betagamma,       "doc"),
        
        
        PROPERTY(nanotron::MCLabel, matchedGenJetDeltaR,  "doc"),
        PROPERTY(nanotron::MCLabel, matchedGenJetPt,      "doc"),
        PROPERTY(nanotron::MCLabel, sharedVertexFraction, "doc"),
        
        PROPERTY(nanotron::MCLabel, genTauMass,  "doc"),
        PROPERTY(nanotron::MCLabel, recoTauMass, "doc")
    };
    
    
    produces<nanoaod::FlatTable>("jetorigin");
}


NANOGenProducer::~NANOGenProducer()
{
}


// ------------ method called to produce the data  ------------
void
NANOGenProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(_jet_src, jets);

    edm::Handle<std::vector<reco::MCLabelInfo>> label_infos;
    iEvent.getByToken(_label_src, label_infos);

    std::size_t njets = jets->size();
    std::size_t ntruth = label_infos->size();


    auto jetOriginTable = std::make_unique<nanoaod::FlatTable>(ntruth, "jetorigin", false, false);



    std::vector<int> truth_jetIdx;

    std::unordered_map<std::string,std::vector<int>> labelList;
    std::unordered_map<std::string,std::vector<int>> tauDecayList;
    

        
    FlatTableFillerList<nanotron::MCLabel> labelPropertyFillerList(labelProperties);
            
    for (std::size_t itag = 0; itag < ntruth; itag++) {
    	const auto& label = label_infos->at(itag).features();
        int jetIdx = -1;
        auto base_jet_ref = label_infos->at(itag).jet();
        if (base_jet_ref.isAvailable() and base_jet_ref.isNonnull()){
        	const auto& base_jet = base_jet_ref.get();
        	for (std::size_t ijet = 0; ijet < njets; ijet++) {
            	auto jet = jets->at(ijet);
            	if (reco::deltaR(base_jet->p4(),jet.p4()) < 1e-4){
                	jetIdx = ijet;
                	break;
            	}
        	}
        } 


        truth_jetIdx.push_back(jetIdx);
        
        for (const auto& labelType: labelTypes) labelList[labelType.first].push_back(label.type == labelType.second ? 1 : 0);
        

        //for (const auto& tauDecayType: tauDecayTypes) tauDecayList[tauDecayType.first].push_back(label.tauDecay == tauDecayType.second ? 1 : 0);
       

        labelPropertyFillerList.push_back(label);
 
    }

    
    jetOriginTable->addColumn<int>("jetIdx", truth_jetIdx, "doc", nanoaod::FlatTable::IntColumn);
    
    for (const auto& labelType: labelTypes) jetOriginTable->addColumn<int>(labelType.first, labelList[labelType.first], "doc", nanoaod::FlatTable::IntColumn);
    
    //for (const auto& tauDecayType: tauDecayTypes) jetOriginTable->addColumn<int>("tauDecay_"+tauDecayType.first, tauDecayList[tauDecayType.first], "doc", nanoaod::FlatTable::IntColumn);

    labelPropertyFillerList.fill(jetOriginTable);
     
    iEvent.put(std::move(jetOriginTable), "jetorigin");
}

void
NANOGenProducer::beginStream(edm::StreamID)
{
}

void
NANOGenProducer::endStream() {
}

 
void
NANOGenProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NANOGenProducer);

