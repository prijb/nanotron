//
//
//
//
//

#ifndef nanotron_DataFormats_MCGhostFlavourInfo_h
#define nanotron_DataFormats_MCGhostFlavourInfo_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "nanotron/DataFormats/interface/MCGenDecayInfo.h"

namespace nanotron {

struct MCGhostFlavour {

    edm::Ptr<nanotron::MCGenDecayInfo> decay;
    edm::Ptr<reco::GenParticle> ghost;
    
    MCGhostFlavour(
        edm::Ptr<nanotron::MCGenDecayInfo> decay,
        edm::Ptr<reco::GenParticle> ghost
    ):
        decay(decay),
        ghost(ghost)
    {
    }
};

struct MCGhostFlavourInfo {
    std::vector<MCGhostFlavour> llpFlavours;
};

}

#endif
