// Neutral object observables
//
//
//
//

#ifndef nanotron_DataFormats_NeutralCandidateFeatures_h
#define nanotron_DataFormats_NeutralCandidateFeatures_h

#include <cmath>

namespace nanotron {

struct NeutralCandidateFeatures  {
    
    float px;
    float py;
    float pz;

    float ptrel;
    float deta;
    float dphi;
    float deltaR;
    
    int isGamma;
    float hcal_fraction;
    float drminsv;
    float puppi_weight;
    float relmassdrop;
    
    NeutralCandidateFeatures():

        px(0),
        py(0),
        pz(0),

        ptrel(0),
        deta(0),
        dphi(0),
        deltaR(0),

        isGamma(0),
        hcal_fraction(0),
        drminsv(0),
        puppi_weight(0),
        relmassdrop(0)
    {}
    
    bool operator<(const NeutralCandidateFeatures& other) const {

        if (std::fabs(drminsv-other.drminsv)>std::numeric_limits<float>::epsilon()) {
            // sort increasing
            return drminsv<other.drminsv;
        }
        else {
            // sort decreasing
            return ptrel>other.ptrel;
        }
        return false;
    }
};

}

#endif
