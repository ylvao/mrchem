/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/Printer>

#include "Functional.h"

namespace mrdft {

Eigen::MatrixXd Functional::evaluate(Eigen::MatrixXd &inp) const {
    int nInp = xc_input_length(*xcfun);  // Input parameters to XCFun
    int nOut = xc_output_length(*xcfun); // Input parameters to XCFun
    int nPts = inp.rows();
    if (nInp != inp.cols()) MSG_ABORT("Invalid input");

    double iDat[nInp];
    double oDat[nOut];

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(nPts, nOut);
    for (int i = 0; i < nPts; i++) {
        bool calc = true;
        for (int j = 0; j < nInp; j++) iDat[j] = inp(i, j);
        if (isSpin()) {
            if (iDat[0] < cutoff and iDat[1] < cutoff) calc = false;
        } else {
            if (iDat[0] < cutoff) calc = false;
        }
        if (calc) {
            xc_eval(*xcfun, iDat, oDat);
            for (int j = 0; j < nOut; j++) out(i, j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i, j) = 0.0;
        }
    }
    return out;
}

Eigen::MatrixXd Functional::contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.rows();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nPts, nFcs);
    out_data.col(0) = xc_data.col(0); // we always keep the energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        Eigen::VectorXd cont_i = Eigen::VectorXd::Zero(nPts);
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            Eigen::VectorXd cont_ij = Eigen::VectorXd::Zero(nPts);
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                cont_ij = xc_data.col(xc_idx).array() * d_data.col(d_idx).array();
            } else {
                cont_ij = xc_data.col(xc_idx);
            }
            cont_i += cont_ij;
        }
        out_data.col(i + 1) = cont_i; // The first column contains the energy functional
    }
    return out_data;
}

} // namespace mrdft
