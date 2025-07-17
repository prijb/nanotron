// Taken from NanoAOD
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <iostream>
#include <algorithm>

class PileupPtHatTableProducer : public edm::global::EDProducer<> {
public:
  PileupPtHatTableProducer(edm::ParameterSet const& params)
  : pileupTag_(consumes<std::vector<PileupSummaryInfo>>(params.getParameter<edm::InputTag>("src"))){
    produces<nanoaod::FlatTable>();
  }

  ~PileupPtHatTableProducer() override {}

  void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
    edm::Handle<std::vector<PileupSummaryInfo>> pileupInfo;
    iEvent.getByToken(pileupTag_, pileupInfo);

    // Loop over pileup info and keep in-time pileup
    std::vector<PileupSummaryInfo> pileupInfoIntime;
    for (const auto& pu: *pileupInfo){
      if (pu.getBunchCrossing() == 0) {
        pileupInfoIntime.push_back(pu);
      }
    }

    std::vector<float> puPtHats;
    if (!pileupInfoIntime.empty()) {
      puPtHats = pileupInfoIntime[0].getPU_pT_hats();
    }
    std::sort(puPtHats.begin(), puPtHats.end(), std::greater<float>());

    // Make the table and fill it
    auto puTable = std::make_unique<nanoaod::FlatTable>(puPtHats.size(), "Pileup", false);
    puTable->addColumn<float>("pthat", puPtHats, "pT hat of pileup events", 20);

    iEvent.put(std::move(puTable)); 
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
  }

  protected:
    const edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupTag_;
};


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PileupPtHatTableProducer);