//
//
//
//
//

#ifndef nanotron_DataFormats_ShallowTagInfoFeatures_h
#define nanotron_DataFormats_ShallowTagInfoFeatures_h

namespace nanotron {

struct ShallowTagInfoFeatures  {

    // Jet general
    float trackSumJetEtRatio;        // Ratio of track sum transverse energy over jet energy
    float trackSumJetDeltaR;         // Pseudoangular distance between jet axis and track fourvector sum
    int vertexCategory;              // Category of secondary vertex (Reco, Pseudo, No)
    float trackSip2dValAboveCharm;   // Track 2D signed impact parameter of first track lifting mass above charm
    float trackSip2dSigAboveCharm;   // Track 2D signed impact parameter significance of first track lifting mass above charm
    float trackSip3dValAboveCharm;   // Track 3D signed impact parameter of first track lifting mass above charm
    float trackSip3dSigAboveCharm;   // Track 3D signed impact parameter significance of first track lifting mass above charm
    
    // Track info
    int jetNTracksEtaRel;            // Tracks associated to jet for which trackEtaRel is calculated
    int jetNSelectedTracks;
    
    ShallowTagInfoFeatures():
        trackSumJetEtRatio(0),       // Ratio of track sum transverse energy over jet energy
        trackSumJetDeltaR(0),        // Pseudoangular distance between jet axis and track fourvector sum
        vertexCategory(0),           // Category of secondary vertex (Reco, Pseudo, No)
        trackSip2dValAboveCharm(0),  // Track 2D signed impact parameter of first track lifting mass above charm
        trackSip2dSigAboveCharm(0),  // Track 2D signed impact parameter significance of first track lifting mass above charm
        trackSip3dValAboveCharm(0),  // Track 3D signed impact parameter of first track lifting mass above charm
        trackSip3dSigAboveCharm(0),  // Track 3D signed impact parameter significance of first track lifting mass above charm
        // track info
        jetNTracksEtaRel(0),         // Tracks associated to jet for which trackEtaRel is calculated
        jetNSelectedTracks(0)
    {}
};

}

#endif
