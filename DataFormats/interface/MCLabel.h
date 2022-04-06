//
//
//
//
//

#ifndef nanotron_DataFormats_MCLabel_h
#define nanotron_DataFormats_MCLabel_h

#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"
#include "nanotron/DataFormats/interface/DisplacedGenVertex.h"

namespace nanotron {

class MCLabel {

    public:
        enum class Type {
            isPU,
            isB,
            isBB,
            isLeptonic_B,
            isC,
            isCC,
            isLeptonic_C,
            isS,
            isUD,
            isG,
            isPrompt_MU,
            isPrompt_E,
            isPrompt_PHOTON,
            isPrompt_TAU,
            
            isMC_RAD,  // no flavour match (likely from wide angle radiation)
            isMC_MU,   // prompt lepton
            isMC_E,    // prompt lepton
            isMC_Q,    // single light quark
            isMC_QMU,  // single light quark + prompt lepton
            isMC_QE,   // single light quark + prompt lepton
            isMC_QQ,   // double light quark
            isMC_QQMU, // double light quark + prompt lepton
            isMC_QQE,  // double light quark + prompt lepton
            isMC_B,    // single b/c quark
            isMC_BMU,  // single b/c quark + prompt lepton
            isMC_BE,   // single b/c quark + prompt lepton
            isMC_BB,   // double b/c quark
            isMC_BBMU, // double b/c quark + prompt lepton
            isMC_BBE,  // double b/c quark + prompt lepton
            isMC_TAU,
            isMC_QTAU,
            isMC_QQTAU,
            isMC_BTAU,
            isMC_BBTAU,
            
            isMC_PHOTON,
            isMC_QPHOTON,
            isMC_QQPHOTON,
            isMC_BPHOTON,
            isMC_BBPHOTON,
            isUndefined
            
        };
        
        Type type;
        
        enum class TauDecay {
            NO_TAU,     // no tau decay
            INVISIBLE,  // tau decay but not reconstructable
            E,          // to electron
            MU,         // to muon
            H,          // 1 charged hadron
            H_1PI0,     // 1 charged hadron + pi0(->2gamma)
            H_XPI0,     // 1 charged hadron + 2 or more pi0(->2gamma)
            HHH,        // 3 charged hadrons
            HHH_XPI0    // 3 charged hadron + 1 or more pi0(->2gamma)
        };
        
        TauDecay tauDecay;
        
        int   partonFlavor;
        int   hadronFlavor;
        int   llpId;
        float llp_mass;
        float llp_pt;
        
        float displacement;
        float displacement_xy;
        float displacement_z;
        float decay_angle;
        float betagamma;
        
        float matchedGenJetDeltaR;
        float matchedGenJetPt;
        float sharedVertexFraction;
        
        float genTauMass;
        float recoTauMass;
        
        MCLabel():
            type(Type::isUndefined),
            tauDecay(TauDecay::NO_TAU),
            partonFlavor(0),
            hadronFlavor(0),
            llpId(0),
            llp_mass(DisplacedGenVertex::MIN_LLP_MASS),
            llp_pt(0),
            
            displacement(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)),    // log10(x / 1cm)
            displacement_xy(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)), // log10(x / 1cm)
            displacement_z(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)),  // log10(x / 1cm)

            decay_angle(0),
            betagamma(0),
            matchedGenJetDeltaR(-1),
            matchedGenJetPt(-1),
            sharedVertexFraction(0),
            genTauMass(-1),
            recoTauMass(-1)
        {
        }
        
        inline static const std::string typeToString(const Type& type) {

            switch (type) {

                case Type::isPU:
                    return "isPU";
                case Type::isB:
                    return "isB";
                case Type::isBB:
                    return "isBB";
                case Type::isLeptonic_B:
                    return "isLeptonic_B";
                case Type::isLeptonic_C:
                    return "isLeptonic_C";
                case Type::isC:
                    return "isC";
                case Type::isCC:
                    return "isCC";
                case Type::isS:
                    return "isS";               
                case Type::isUD:
                    return "isUD";
                case Type::isG:
                    return "isG";

                case Type::isPrompt_MU:
                    return "isPrompt_MU";
                case Type::isPrompt_E:
                    return "isPrompt_E";
                case Type::isPrompt_TAU:
                    return "isPrompt_TAU";
                case Type::isPrompt_PHOTON:
                    return "isPrompt_PHOTON";
                
                case Type::isMC_RAD:
                    return "isMC_RAD";
                    
                case Type::isMC_MU:
                    return "isMC_MU";
                case Type::isMC_E:
                    return "isMC_E";
                    
                case Type::isMC_Q:
                    return "isMC_Q";
                    
                case Type::isMC_QMU:
                    return "isMC_QMU";
                case Type::isMC_QE:
                    return "isMC_QE";
                    
                case Type::isMC_QQ:
                    return "isMC_QQ";
                    
                case Type::isMC_QQMU:
                    return "isMC_QQMU";
                case Type::isMC_QQE:
                    return "isMC_QQE"; 
                    
                case Type::isMC_B:
                    return "isMC_B";
                    
                case Type::isMC_BMU:
                    return "isMC_BMU";
                case Type::isMC_BE:
                    return "isMC_BE";
                    
                case Type::isMC_BB:
                    return "isMC_BB";
                    
                case Type::isMC_BBMU:
                    return "isMC_BBMU";
                case Type::isMC_BBE:
                    return "isMC_BBE";

                case Type::isMC_TAU:
                    return "isMC_TAU";
                case Type::isMC_QTAU:
                    return "isMC_QTAU";
                case Type::isMC_QQTAU:
                    return "isMC_QQTAU";
                case Type::isMC_BTAU:
                    return "isMC_BTAU";
                case Type::isMC_BBTAU:
                    return "isMC_BBTAU";
                    
                case Type::isMC_PHOTON:
                    return "isMC_TAU";
                case Type::isMC_QPHOTON:
                    return "isMC_QTAU";
                case Type::isMC_QQPHOTON:
                    return "isMC_QQTAU";
                case Type::isMC_BPHOTON:
                    return "isMC_BTAU";
                case Type::isMC_BBPHOTON:
                    return "isMC_BBTAU";
                    
                case Type::isUndefined:
                    return "isUndefined";
            }
            return "isUndefined";
        }
};

}

#endif
