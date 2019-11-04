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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/trees/FunctionNode.h>

#include "xc_utils.h"

namespace mrdft {

namespace xc_utils {
void fill_output_mask(Eigen::MatrixXi &mask, int start);
} // namespace xc_utils

Eigen::MatrixXi xc_utils::build_output_mask(bool is_lda, bool is_spin_sep, int order) {
    int start = 2;
    bool is_gga = not is_lda;
    Eigen::MatrixXi mask(1, 1);
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

Eigen::VectorXi xc_utils::build_density_mask(bool is_lda, bool is_spin_sep, int order) {
    bool is_gga = not is_lda;
    Eigen::VectorXi mask(1);
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

void xc_utils::fill_output_mask(Eigen::MatrixXi &mask, int value) {
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

std::vector<mrcpp::FunctionNode<3> *> xc_utils::fetch_nodes(int n, mrcpp::FunctionTreeVector<3> &inp) {
    std::vector<mrcpp::FunctionNode<3> *> nodes;
    for (auto i = 0; i < inp.size(); i++) {
        auto &iTree = mrcpp::get_func(inp, i);
        auto &iNode = iTree.getEndFuncNode(n);
        nodes.push_back(&iNode);
    }
    return nodes;
}

/** @brief Converts data from a FunctionNode to a matrix
 *
 * The FunctionNode(s) row data is packed into a matrix whose
 * dimensions are the overall number of grid points (nCoefs) and the
 * number of functions (nFuncs).
 *
 * parma[in] n the Index of the requested node
 * param[in] nFuncs The number of functions
 * param[in] trees The array of FunctionTree(s)
 * param[in] data The matrix object
 */
Eigen::MatrixXd xc_utils::compress_nodes(std::vector<mrcpp::FunctionNode<3> *> &inp_nodes) {
    Eigen::MatrixXd out_data;
    auto nFuncs = inp_nodes.size();
    if (nFuncs > 0) {
        auto nCoefs = inp_nodes[0]->getNCoefs();
        out_data = Eigen::MatrixXd::Zero(nCoefs, nFuncs);
        for (auto i = 0; i < nFuncs; i++) {
            auto &node = inp_nodes[i];
            Eigen::VectorXd col_i;
            node->getValues(col_i);
            if (col_i.size() != nCoefs) MSG_ABORT("Size mismatch");
            out_data.col(i) = col_i;
        }
    }
    return out_data;
}

/** @brief Converts data from a matrix to a FunctionNode
 *
 * The matrix containing the output from xcfun is converted back to the corresponding FunctionNode(s). The matrix
 * dimensions are the overall number of grid points (nCoefs) and the number of functions (nFuncs).
 *
 * parma[in] n the Index of the requested node
 * param[in] nFuncs The number of functions
 * param[in] trees The array of FunctionTree(s)
 * param[in] data The matrix object
 */
void xc_utils::expand_nodes(std::vector<mrcpp::FunctionNode<3> *> &out_nodes, Eigen::MatrixXd &out_data) {
    auto nFuncs = out_nodes.size();
    if (out_data.cols() != nFuncs) MSG_ERROR("Size mismatch");

    for (auto i = 0; i < nFuncs; i++) {
        auto &node = out_nodes[i];
        node->setValues(out_data.col(i));
    }
}

mrcpp::FunctionTreeVector<3> xc_utils::log_gradient(mrcpp::DerivativeOperator<3> &diff_oper,
                                                    mrcpp::FunctionTree<3> &rho) {
    mrcpp::FunctionTree<3> zeta(rho.getMRA());
    mrcpp::copy_grid(zeta, rho);
    mrcpp::copy_func(zeta, rho);
    for (auto i = 0; i < zeta.getNEndNodes(); i++) {
        auto &node = zeta.getEndFuncNode(i);
        Eigen::VectorXd values;
        node.getValues(values);
        for (auto j = 0; j < node.getNCoefs(); j++) {
            if (values[j] > mrcpp::MachineZero) {
                values[j] = std::log(values[j]);
            } else {
                values[j] = mrcpp::MachineZero;
            }
        }
        node.setValues(values);
    }
    zeta.mwTransform(mrcpp::BottomUp);

    mrcpp::FunctionTreeVector<3> grad_zeta = mrcpp::gradient(diff_oper, zeta);

    mrcpp::FunctionTreeVector<3> grad_rho;
    for (int i = 0; i < 3; i++) {
        mrcpp::FunctionTree<3> *grad_comp = new mrcpp::FunctionTree<3>(rho.getMRA());
        mrcpp::copy_grid(*grad_comp, rho);
        mrcpp::multiply(-1.0, *grad_comp, 1.0, rho, mrcpp::get_func(grad_zeta, i));
        grad_rho.push_back(std::make_tuple(1.0, grad_comp));
    }
    mrcpp::clear(grad_zeta, true);
    return grad_rho;
}

} // namespace mrdft
