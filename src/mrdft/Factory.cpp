#include "Factory.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include "LibXC.h"
// #include "LibXC.cpp"


#include "GGA.h"
#include "Grid.h"
#include "LDA.h"
#include "MRDFT.h"
#include "SpinGGA.h"
#include "SpinLDA.h"

namespace mrdft {

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
        : mra(MRA) {}

// Factory::~Factory() {
//     cleanupFunctionals();
// }

void Factory::cleanupFunctionals() {
    for (auto &func_data : functionals) {
        if (func_data.initialized) {
            xc_func_end(&func_data.func);
            func_data.initialized = false;
        }
    }
    functionals.clear();
}

// int Factory::mapFunctionalName(const std::string &name) const {
//     // Map common functional names to LibXC IDs
//     if (name == "LDA" || name == "LDA_X") return XC_LDA_X;
//     if (name == "VWN" || name == "LDA_C_VWN") return XC_LDA_C_VWN;
//     if (name == "PBE_X") return XC_GGA_X_PBE;
//     if (name == "PBE_C") return XC_GGA_C_PBE;
//     if (name == "B88") return XC_GGA_X_B88;
//     if (name == "LYP") return XC_GGA_C_LYP;
//     if (name == "B3LYP") return XC_HYB_GGA_XC_B3LYP;
    
//     // If not a common name, try to get it directly from LibXC
//     int func_id = xc_functional_get_number(name.c_str());
//     if (func_id <= 0) {
//         std::string msg = "Unknown functional: " + name;
//         MSG_ABORT(msg.c_str());
//     }
//     return func_id;
// }

// void Factory::setFunctional(const std::string &name, double weight) {
//     int func_id = mapFunctionalName(name);
    
//     LibXCData func_data;
//     func_data.func_id = func_id;
//     func_data.weight = weight;
//     func_data.initialized = false;
    
//     functionals.push_back(func_data);
// }

bool Factory::isGGA() const {
    for (const auto &func_data : functionals) {
        const xc_func_type *func_ptr = &func_data.func;
        if (func_data.initialized) {
            if (func_ptr->info->family == XC_FAMILY_GGA || 
                func_ptr->info->family == XC_FAMILY_HYB_GGA) {
                return true;
            }
        } else {
            // If not initialized, check the family based on ID
            xc_func_type temp_func;
            if (xc_func_init(&temp_func, func_data.func_id, XC_UNPOLARIZED) == 0) {
                bool is_gga = (temp_func.info->family == XC_FAMILY_GGA || 
                              temp_func.info->family == XC_FAMILY_HYB_GGA);
                xc_func_end(&temp_func);
                if (is_gga) return true;
            }
        }
    }
    return false;
}

// bool Factory::isHybrid() const {
//     for (const auto &func_data : functionals) {
//         if (func_data.initialized) {
//             if (xc_hyb_type(&func_data.func) != XC_HYB_NONE) {
//                 return true;
//             }
//         } else {
//             // Check if it's a hybrid by temporarily initializing
//             xc_func_type temp_func;
//             if (xc_func_init(&temp_func, func_data.func_id, XC_UNPOLARIZED) == 0) {
//                 bool is_hybrid = (xc_hyb_type(&temp_func) != XC_HYB_NONE);
//                 xc_func_end(&temp_func);
//                 if (is_hybrid) return true;
//             }
//         }
//     }
//     return false;
// }

bool Factory::isHybrid() const {
    for (const auto &func_data : functionals) {
        const xc_func_type *func_ptr = &func_data.func;
        if (func_data.initialized) {
            if (func_ptr->info->family == XC_FAMILY_HYB_GGA ||
                func_ptr->info->family == XC_FAMILY_HYB_MGGA) {
                return true;
            }
        } else {
            xc_func_type temp_func;
            if (xc_func_init(&temp_func, func_data.func_id, XC_UNPOLARIZED) == 0) {
                bool is_hybrid = (temp_func.info->family == XC_FAMILY_HYB_GGA ||
                                  temp_func.info->family == XC_FAMILY_HYB_MGGA);
                xc_func_end(&temp_func);
                if (is_hybrid) return true;
            }
        }
    }
    return false;
}



// double Factory::getHybridCoeff() const {
//     double coeff = 0.0;
//     for (const auto &func_data : functionals) {
//         if (func_data.initialized) {
//             if (xc_hyb_type(&func_data.func) != XC_HYB_NONE) {
//                 coeff += func_data.weight * xc_hyb_exx_coef(&func_data.func);
//             }
//         } else {
//             // Get hybrid coefficient by temporarily initializing
//             xc_func_type temp_func;
//             if (xc_func_init(&temp_func, func_data.func_id, XC_UNPOLARIZED) == 0) {
//                 if (xc_hyb_type(&temp_func) != XC_HYB_NONE) {
//                     coeff += func_data.weight * xc_hyb_exx_coef(&temp_func);
//                 }
//                 xc_func_end(&temp_func);
//             }
//         }
//     }
//     return coeff;
// }


double Factory::getHybridCoeff() const {
    double coeff = 0.0;
    for (const auto &func_data : functionals) {
        if (func_data.initialized) {
            if (func_data.func.info->family == XC_FAMILY_HYB_GGA ||
                func_data.func.info->family == XC_FAMILY_HYB_MGGA) {
                coeff += func_data.weight * xc_hyb_exx_coef(&func_data.func);
            }
        } else {
            xc_func_type temp_func;
            if (xc_func_init(&temp_func, func_data.func_id, XC_UNPOLARIZED) == 0) {
                if (temp_func.info->family == XC_FAMILY_HYB_GGA ||
                    temp_func.info->family == XC_FAMILY_HYB_MGGA) {
                    coeff += func_data.weight * xc_hyb_exx_coef(&temp_func);
                }
                xc_func_end(&temp_func);
            }
        }
    }
    return coeff;
}


/** @brief Build a MRDFT object from the currently defined parameters */
std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);
    
    // Initialize LibXC functionals
    for (auto &func_data : functionals) {
        int polarization = spin ? XC_POLARIZED : XC_UNPOLARIZED;
        if (xc_func_init(&func_data.func, func_data.func_id, polarization) != 0) {
            std::string msg = "Error initializing LibXC functional";
            MSG_ABORT(msg.c_str());
        }
        func_data.initialized = true;
    }
    
    // Check if we have any functionals
    if (functionals.empty()) {
        MSG_ABORT("No functionals defined");
    }
    
    // Init MW derivative
    bool gga = isGGA();
    if (gga) {
        if (diff_s == "bspline") diff_p = std::make_unique<mrcpp::BSOperator<3>>(mra, 1);
        if (diff_s == "abgv_00") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        if (diff_s == "abgv_55") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
    }
    
    // Had some error with the way these constructors are defined -- Dont need spin yet
    
    // Init XC functional
    std::unique_ptr<Functional> func_p{nullptr};
    bool lda = !gga;
    // if (spin) {
    //     if (gga) func_p = std::make_unique<SpinGGA>(order, functionals, diff_p);
    //     if (lda) func_p = std::make_unique<SpinLDA>(order, functionals);
    // } else {
    //     if (gga) func_p = std::make_unique<GGA>(order, functionals, diff_p);
    //     if (lda) func_p = std::make_unique<LDA>(order, functionals);
    // }
    
    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);
    
    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft