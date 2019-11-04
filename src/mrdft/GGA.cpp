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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"

#include "GGA.h"
#include "xc_utils.h"

namespace mrdft {

GGA::GGA(int k, std::unique_ptr<xc_functional> &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d, bool lg)
        : Functional(k, f)
        , log_grad(lg)
        , derivative(std::move(d)) {
    xc_mask = xc_utils::build_output_mask(false, false, this->order);
    d_mask = xc_utils::build_density_mask(false, false, this->order);
}

void GGA::clear() {
    mrcpp::clear(this->rho, false);
    mrcpp::clear(this->grad, true);
}

int GGA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 4;
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 */
mrcpp::FunctionTreeVector<3> GGA::setupXCInput() {
    if (this->rho.size() < 1) MSG_ERROR("Density not initialized");
    if (this->grad.size() < 3) MSG_ERROR("Gradient not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho[0]);
    out_vec.insert(out_vec.end(), this->grad.begin(), this->grad.begin() + 3);
    return out_vec;
}

/** @brief Allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 */
mrcpp::FunctionTreeVector<3> GGA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (this->order == 2) {
        out_vec.push_back(this->rho[1]);
        out_vec.insert(out_vec.end(), this->grad.begin() + 3, this->grad.begin() + 6);
    }
    return out_vec;
}

/** @brief Allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 */
void GGA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != this->order) MSG_ERROR("Invalid input length");
    if (this->rho.size() > 0) MSG_ERROR("Density not empty");
    if (this->grad.size() > 0) MSG_ERROR("Gradient not empty");

    int n = 0;
    for (int i = 0; i < this->order; i++) this->rho.push_back(inp_vec[n++]);

    for (int i = 0; i < this->order; i++) {
        mrcpp::FunctionTreeVector<3> tmp;
        if (this->log_grad and i == 0) {
            tmp = xc_utils::log_gradient(*this->derivative, mrcpp::get_func(this->rho, i));
        } else {
            tmp = mrcpp::gradient(*this->derivative, mrcpp::get_func(this->rho, i));
        }
        this->grad.insert(this->grad.end(), tmp.begin(), tmp.end());
    }
}

/** @brief Potential calculation for GGA functionals
 *
 * The potential functions are assembled from the xcfun output functions.
 * The method used depends on whether the functional is spin-separated
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 */
mrcpp::FunctionTreeVector<3> GGA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Energy density
    mrcpp::FunctionTree<3> &f_xc = mrcpp::get_func(inp_vec, 0);
    inp_vec[0] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Potential part
    mrcpp::FunctionTree<3> &df_dr = mrcpp::get_func(inp_vec, 1);
    mrcpp::FunctionTreeVector<3> df_dg(inp_vec.begin() + 2, inp_vec.begin() + 5);

    auto *tmp = new mrcpp::FunctionTree<3>(df_dr.getMRA());
    mrcpp::divergence(*tmp, *this->derivative, df_dg);

    auto *v_xc = new mrcpp::FunctionTree<3>(df_dr.getMRA());
    mrcpp::build_grid(*v_xc, df_dr);
    mrcpp::build_grid(*v_xc, *tmp);
    mrcpp::add(-1.0, *v_xc, 1.0, df_dr, -1.0, *tmp);
    delete tmp;

    // Collect output
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, v_xc));
    v_xc = nullptr;

    return out_vec;
}

} // namespace mrdft
