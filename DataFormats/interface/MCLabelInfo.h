//
//
//
//
//

#ifndef nanotron_DataFormats_MCLabelInfo_h
#define nanotron_DataFormats_MCLabelInfo_h

#include "nanotron/DataFormats/interface/MCLabel.h"
#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

namespace reco {

typedef FeaturesTagInfo<nanotron::MCLabel> MCLabelInfo;

DECLARE_EDM_REFS(MCLabelInfo)

}

#endif
