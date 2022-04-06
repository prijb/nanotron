//
//
//
//
//

#ifndef nanotron_DataFormats_OnionFeatures_h
#define nanotron_DataFormats_OnionFeatures_h


#include <vector>

#include "nanotron/DataFormats/interface/JetFeatures.h"
#include "nanotron/DataFormats/interface/SecondaryVertexFeatures.h"
#include "nanotron/DataFormats/interface/ShallowTagInfoFeatures.h"
#include "nanotron/DataFormats/interface/NeutralCandidateFeatures.h"
#include "nanotron/DataFormats/interface/ChargedCandidateFeatures.h"
#include "nanotron/DataFormats/interface/MuonCandidateFeatures.h"
#include "nanotron/DataFormats/interface/ElectronCandidateFeatures.h"


namespace nanotron {

class OnionTagFeatures {
    
  public:
    JetFeatures jet_features;
    
    ShallowTagInfoFeatures tag_info_features;
    std::vector<SecondaryVertexFeatures> sv_features;
    std::vector<SecondaryVertexFeatures> sv_adapted_features;

    std::vector<NeutralCandidateFeatures> npf_features;
    std::vector<ChargedCandidateFeatures> cpf_features;

    std::vector<MuonCandidateFeatures>     mu_features;
    std::vector<ElectronCandidateFeatures> elec_features;
    
    std::size_t npv;
};

}  

#endif
