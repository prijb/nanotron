//
//
//
//
//

#ifndef nanotron_DataFormats_OnionInfo_h
#define nanotron_DataFormats_OnionInfo_h

#include "nanotron/DataFormats/interface/OnionTagFeatures.h"
#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

namespace reco {

    typedef FeaturesTagInfo<nanotron::OnionTagFeatures> OnionTagInfo;
    
    DECLARE_EDM_REFS(OnionTagInfo)
}

#endif

