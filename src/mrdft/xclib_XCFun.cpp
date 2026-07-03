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
#include "xclib_XCFun.h"
#include "xc_func_alias.h"

namespace mrdft {

double XCFun::setFunctional(const std::string &name, double c, double cutoff, bool spin) {
    double customExx = 0.0;
    xcfun_set(xcfun, name.c_str(), c);
    return customExx;
}

void XCFun::initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, bool spin, int order, bool gamma) {
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

double XCFun::getAmountExx() const {
    double exx = 0;
    xcfun_get(xcfun, "exx", &exx);
    return exx;
}

void XCFun::printFunctionalReference(int out_txt_width, std::vector<std::string> xcfun_func_names) const {
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
        printout(0, xcfun_splash());
        std::cout << "\nXCFun functionals used in this calculation:\n";
        for (const auto &func_name : xcfun_func_names) {
            std::string xcfun_ref = xcfun_describe_long(func_name.c_str());
            std::string xcfun_ref_str = "  - " + xcfun_ref;
            print_wrap(xcfun_ref_str, out_txt_width, 4);
        return;
    }
}

void XCFun::callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, bool spin, double cutoff) const {
    if (nInp != xcfun_input_length(xcfun) or nOut != xcfun_output_length(xcfun)) { throw std::logic_error("Dimension mismatch!\n"); }

    for (int i = 0; i < nPts; i++) {
        if (spin) {
            if (inp(0, i) < cutoff and inp(1, i) < cutoff) continue;
        } else {
            if (inp(0, i) < cutoff) continue;
        }
        xcfun_eval(xcfun, inp.col(i).data(), out.col(i).data());
    }
}

} // mrdft
