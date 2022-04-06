//
//
//
//
//

#ifndef nanotron_DataFormats_MCGenDecayInfo_h
#define nanotron_DataFormats_MCGenDecayInfo_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

namespace nanotron {

struct MCGenDecayInfo {
    std::string name;
    edm::Ptr<reco::GenParticle> llp;
    edm::PtrVector<reco::GenParticle> decayProducts;
};



}

#endif
