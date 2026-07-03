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

#include <memory>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>
#include <xc_funcs.h>
#include <xc.h>
#include "xclib.h"
#include "xc_func_alias.h"

namespace mrdft {

double XClib::setFunctional(const std::string &name, double c, double cutoff, bool spin) {
    double customExx = 0.0;
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
        xcfun_set(xcfun, name.c_str(), c);
    }
    return customExx;
}

void XClib::initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, bool spin, int order, bool gamma) {
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
                    MSG_ABORT("Meta-GGA functionals are not supported in MRChem.!\n");
                    mgga = true;

                default:
                    MSG_ABORT("Libxc functional family not handled in MRChem!\n");
            }

            // Check if functional is range separated
            if ((f->info->flags & XC_FLAGS_HYB_CAM)  || (f->info->flags & XC_FLAGS_HYB_LC)) MSG_ABORT("Coulomb attenuated functionals not supported in MRChem!\n");
            if ((f->info->flags & XC_FLAGS_HYB_CAMY) || (f->info->flags & XC_FLAGS_HYB_LCY)) MSG_ABORT("Yukawa attenuated functionals not supported in MRChem!\n");
        }
    } else {
        gga = xcfun_is_gga(xcfun);
        lda = not gga;
        unsigned int mode = 1;                    //!< only partial derivative mode implemented
        unsigned int func_type = (gga) ? 1 : 0;   //!< only LDA and GGA supported for now
        unsigned int dens_type = 1 + spin;        //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
        unsigned int laplacian = 0;               //!< no laplacian
        unsigned int kinetic = 0;                 //!< no kinetic energy density
        unsigned int current = 0;                 //!< no current density
        unsigned int exp_derivative = not(gamma); //!< use gamma or explicit derivatives
        if (not(gga)) exp_derivative = 0;         //!< fall back to gamma-type derivatives if LDA
        xcfun_user_eval_setup(xcfun, order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);
    }
}

double XClib::getAmountExx() const {
    double exx = 0;
    if (libxc) {
        for (std::size_t i = 0; i < libxc_objects.size(); ++i) {
            const xc_func_type *f = libxc_objects[i];
            double frac = xc_hyb_exx_coef(f);
            exx += libxc_coefs[i] * frac;
        }
    } else {
        xcfun_get(xcfun, "exx", &exx);
    }
    return exx;
}

void XClib::printFunctionalReference(int out_txt_width, std::vector<std::string> xcfun_func_names) const {
    // Only run on main thread
    if (mrcpp::mpi::world_rank != 0) {
        return;
    }

    auto pwidth = mrcpp::Printer::getWidth();
    int txt_width = 50;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << "                  XC Functional                   " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    mrcpp::print::separator(0, '*', 1);


    // Conditional reference printing
    auto print_wrap = [&](std::string str, std::size_t txt_width, int indent = 0) {
        const std::string continuation_indent(indent, ' ');
        size_t offset = 0;
        while (offset + txt_width < str.size()) {
            // Is the line already sufficiently short?
            size_t newline_pos = str.find('\n', offset);
            if (newline_pos < offset + txt_width) {
                if (newline_pos != std::string::npos && newline_pos + 1 < str.size()) {
                    str.insert(newline_pos + 1, continuation_indent);
                    offset = newline_pos + 1 + continuation_indent.size();
                } else {
                    offset = newline_pos + 1;
                }
                continue;
            }

            size_t space_pos = str.rfind(' ', offset + txt_width);
            // If the string doesn't have a space, or it is too far out, hard insert a newline
            if (space_pos == std::string::npos || space_pos - offset > txt_width) {
                space_pos = offset + txt_width;
                str.insert(space_pos, "\n" + continuation_indent);
            } else {
                str[space_pos] = '\n';
                str.insert(space_pos + 1, continuation_indent);
            }
            offset = space_pos + 1 + continuation_indent.size();
        }
        std::cout << str;
    };

    if (libxc) {
        std::string libxc_ref_str = "Using Libxc (version " + std::string(xc_version_string()) + ") to evaluate density functionals. Libxc is free software. It is " +
                                    "distributed under the Mozilla Public License, version 2.0. For " + "more information, please check the Libxc manual. You should cite\n\n" +
                                    xc_reference() + " DOI: " + xc_reference_doi() + "\n\nwhen " + "reporting the results of your calculation in a scientific article.\n";
        print_wrap(libxc_ref_str, out_txt_width);

        // Avoid printing the same LibXC functional multiple times
        std::set<int> printed_ids;
        std::cout << "\nLibxc functionals used in this calculation:\n";
        for (const auto *func : libxc_objects) {
            if (func->info == nullptr) continue; // safety

            int id = func->info->number;
            if (!printed_ids.insert(id).second) {
                // Already printed this ID
                continue;
            }

            char *name = xc_functional_get_name(id);
            std::string func_id_str = "  - " + std::string(name) + " (ID " + std::to_string(id) + "): " + func->info->name + "\n";
            print_wrap(func_id_str, out_txt_width);
            free(name);

            for (int number = 0; number < XC_MAX_REFERENCES; number++) {
                auto reference = xc_func_info_get_references(func->info, number);
                if (reference == nullptr) break;
                std::string func_ref_str = "     * " + std::string(xc_func_reference_get_ref(reference)) + ", DOI:" + xc_func_reference_get_doi(reference) + "\n";
                print_wrap(func_ref_str, out_txt_width, 7);
            }
        }
    }
    if (not XClib::libxc) {
        printout(0, xcfun_splash());
        std::cout << "\nXCFun functionals used in this calculation:\n";
        for (const auto &func_name : xcfun_func_names) {
            std::string xcfun_ref = xcfun_describe_long(func_name.c_str());
            std::string xcfun_ref_str = "  - " + xcfun_ref;
            print_wrap(xcfun_ref_str, out_txt_width, 4);
        }
        return;
    }
}

}; // mrdft
