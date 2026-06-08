/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "Factory.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>

#include "MRDFT.h"
#include "Grid.h"
#include "LDA.h"
#include "GGA.h"
#include "mGGA.h"
#include "SpinGGA.h"
#include "SpinLDA.h"
#include "SpinmGGA.h"

namespace mrdft {

bool Factory::libxc = false;

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
        : mra(MRA)
        , xcfun_p(xcfun_new(), xcfun_delete) {}

void Factory::setFunctional(const std::string &name, double c) {
    setLibxc(libxc); // should probably be where setFunctional is called

    if (libxc) {
        std::vector<int> ids;
        std::vector<double> coefs;
        double scaled_custom_exx = 0.0;

        mapFunctionalName(name, ids, coefs, scaled_custom_exx);
        customExx += c * scaled_custom_exx;
        xc_func_type *libxc_obj;
        for (size_t i = 0; i < ids.size(); i++) {
            libxc_obj = xc_func_alloc();
            auto return_code = xc_func_init(libxc_obj, ids[i], spin ? XC_POLARIZED : XC_UNPOLARIZED);
            if (return_code != 0) {
                std::ostringstream oss;
                oss << "Functional id = " << ids[i] << " not found in the employed version of Libxc.\n";
                MSG_ABORT(oss.str());
            }
            xc_func_set_dens_threshold(libxc_obj, cutoff);
            libxc_objects.push_back(libxc_obj);
            libxc_coefs.push_back(c * coefs[i]);
        }

    } else {
        xcfun_set(xcfun_p.get(), name.c_str(), c);
        xcfun_func_names.push_back(name);
    }
}

std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);
    setLibxc(libxc);

    // Init XCFun or Libxc
    bool lda = false;
    bool gga = false;
    bool mgga = false;
    if (libxc) {
        for (const auto *f : libxc_objects) {

            switch (f->info->family) {
                case XC_FAMILY_LDA:
#ifdef XC_FAMILY_HYB_GGA
                case XC_FAMILY_HYB_LDA:
#endif
                    lda = true;
                    break;

                case XC_FAMILY_GGA:
#ifdef XC_FAMILY_HYB_GGA
                case XC_FAMILY_HYB_GGA:
#endif
                    gga = true;
                    break;

                case XC_FAMILY_MGGA:
                case XC_FAMILY_HYB_MGGA:
                    mgga = true;
                    break;

                default:
                    MSG_ABORT("Libxc functional family not handled in MRChem!\n");
            }

            // Check if functional is range separated
            if ((f->info->flags & XC_FLAGS_HYB_CAM)  || (f->info->flags & XC_FLAGS_HYB_LC)) MSG_ABORT("Coulomb attenuated functionals not supported in MRChem!\n");
            if ((f->info->flags & XC_FLAGS_HYB_CAMY) || (f->info->flags & XC_FLAGS_HYB_LCY)) MSG_ABORT("Yukawa attenuated functionals not supported in MRChem!\n");
        }
    } else {
        gga = xcfun_is_gga(xcfun_p.get());
        mgga = xcfun_is_metagga(xcfun_p.get());
        if (not(gga)) { if (not(mgga)) {lda = true; } }
        unsigned int mode = 1;                    //!< only partial derivative mode implemented
        unsigned int dens_type = 1 + spin;        //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
        unsigned int laplacian = 0;               //!< no laplacian
        unsigned int kinetic = mgga ? 1u : 0u;    //!< request kinetic energy density for meta-GGAs
        unsigned int current = 0;                 //!< no current density
        unsigned int exp_derivative = not(gamma); //!< use gamma or explicit derivatives
        if (lda) exp_derivative = 0;              //!< fall back to gamma-type derivatives if LDA
        unsigned int func_type = 0;               //!< LDA = 0, GGA = 1, mGGA = 2
        if (mgga) {
            func_type = 2;
        } else if (gga) {
            func_type = 1;
        } else {
            func_type = 0;
        }
        xcfun_user_eval_setup(xcfun_p.get(), order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);
    }

    // bool lda = not gga;

    // Init MW derivative
    if (gga) {
        // make shared pointer???
        if (diff_s == "bspline") diff_p = std::make_shared<mrcpp::BSOperator<3>>(mra, 1);
        if (diff_s == "abgv_00") diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        if (diff_s == "abgv_55") diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
    }

    if (mgga) {
        kin_p = std::make_unique<mrchem::KineticOperator>(mrchem::KineticOperator(diff_p));
    }
    std::unique_ptr<Functional> func_p{nullptr};
    if (spin) {
        if (mgga) {
            // TODO: Spin-mGGA class analogous to SpinGGA
            // MSG_ABORT("Spin-polarized meta-GGA class not implemented yet.");
            func_p = std::make_unique<SpinmGGA>(order, xcfun_p, diff_p);
        } else if (gga) {
            func_p = std::make_unique<SpinGGA>(order, xcfun_p, diff_p);
        } else if (lda) {
            func_p = std::make_unique<SpinLDA>(order, xcfun_p);
        } else {
            MSG_ABORT("Case not handled");
        }
    } else {
        if (mgga) {
            func_p = std::make_unique<mGGA>(order, xcfun_p, diff_p);
        } else if (gga) {
            func_p = std::make_unique<GGA>(order, xcfun_p, diff_p);
        } else if (lda) {
            func_p = std::make_unique<LDA>(order, xcfun_p);
        } else {
            MSG_ABORT("Case not handled");
        }
    }

    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    if (libxc) { func_p->setLibxcFunctionalObject(libxc_objects, libxc_coefs); }
    if (not libxc) { func_p->setXCFunFunctionalNames(xcfun_func_names); }
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setCustomExx(customExx);
    // Is this ever used?????:
    func_p->setDerivOp(diff_p);
    func_p->setKinOp(kin_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
