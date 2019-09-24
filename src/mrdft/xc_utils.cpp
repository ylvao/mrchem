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

#include "MRCPP/Printer"

#include "xc_utils.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

static int xc_iteration = 0;

using namespace mrdft;

MatrixXi xc_utils::build_output_mask(bool is_lda, bool is_spin_sep, int order) {
    int start = 2;
    bool is_gga = not is_lda;
    MatrixXi mask(1, 1);
    mask << 1;
    switch (order) {
        case 0:
            break;
        case 1:
            if (is_lda and is_spin_sep) {
                mask.resize(2, 1);
                mask << 1, 2;
            } else if (is_gga and not is_spin_sep) {
                mask.resize(4, 1);
                mask << 1, 2, 3, 4;
            } else if (is_gga and is_spin_sep) {
                mask.resize(8, 1);
                mask << 1, 2, 3, 4, 5, 6, 7, 8;
            }
            break;
        case 2:
            if (is_lda and is_spin_sep) {
                start = 3;
                mask.resize(2, 2);
            } else if (is_gga and not is_spin_sep) {
                start = 5;
                mask.resize(4, 4);
            } else if (is_gga and is_spin_sep) {
                start = 9;
                mask.resize(8, 8);
            }
            fill_output_mask(mask, start);
            break;
        default:
            MSG_ABORT("Not implemented");
    }
    return mask;
}

VectorXi xc_utils::build_density_mask(bool is_lda, bool is_spin_sep, int order) {
    bool is_gga = not is_lda;
    VectorXi mask(1);
    switch (order) {
        case 0:
        case 1:
            mask(0) = -1;
            break;
        case 2:
            mask(0) = 0;
            if (is_lda and is_spin_sep) {
                mask.resize(2);
                mask << 0, 1;
            } else if (is_gga and not is_spin_sep) {
                mask.resize(4);
                mask << 0, 1, 2, 3;
            } else if (is_gga and is_spin_sep) {
                mask.resize(8);
                mask << 0, 1, 2, 3, 4, 5, 6, 7;
            }
            break;
        default:
            MSG_ABORT("Not implemented");
    }
    return mask;
}

void xc_utils::fill_output_mask(MatrixXi &mask, int value) {
    for (int i = 0; i < mask.rows(); i++) {
        mask(i, i) = value;
        value++;
        for (int j = i + 1; j < mask.cols(); j++) {
            mask(i, j) = value;
            mask(j, i) = value;
            value++;
        }
    }
}
