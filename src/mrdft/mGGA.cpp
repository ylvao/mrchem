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
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"

#include "mGGA.h"
#include "xc_utils.h"

namespace mrdft {

mGGA::mGGA(int k, XC_p &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> &d)
        : Functional(k, f)
        , derivative(d) {
    xc_mask = xc_utils::build_output_mask(false, false, this->order);
    d_mask = xc_utils::build_density_mask(false, false, this->order);
}

/** @brief Clear internal functions
 *
 * Ownership of densities is outside MRDFT -> clear
 * Ownership of gradients and tau are inside MRDFT -> free
 */
void mGGA::clear() {
    mrcpp::clear(this->rho, false);
    mrcpp::clear(this->grad, true);
    mrcpp::clear(this->tau, true);
}

/** @brief Number of function involved in contraction step */
int mGGA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 5;  // rho_1 + 3 * grad(rho_1) + tau_1
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Collect input functions to xcfun / Libxc evaluation step
 *
 * For mGGA we provide the density part with tau:
 *   [rho_0, grad(rho_0), tau_0]
 *
 * The orbital kinetic energy density tau is computed independently
 * (from the orbitals in XCPotentialD1) and included in the input.
 */
mrcpp::FunctionTreeVector<3> mGGA::setupXCInput() {
    if (this->rho.size() < 1) MSG_ERROR("Density not initialized");
    if (this->grad.size() < 3) MSG_ERROR("Gradient not initialized");
    if (this->tau.size() < 1) MSG_ERROR("Kinetic energy density not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho[0]);
    out_vec.insert(out_vec.end(), this->grad.begin(), this->grad.begin() + 3);
    out_vec.push_back(this->tau[0]);
    return out_vec;
}

/** @brief Collect input functions to contraction step
 *
 * For mGGA:
 * Ground State: No contraction, empty vector
 * Linear Response: [rho_1, grad(rho_1), tau_1]
 * Higher Response: NOT_IMPLEMENTED
 */
mrcpp::FunctionTreeVector<3> mGGA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (this->order == 2) {
        if (this->rho.size() < 2) MSG_ERROR("Perturbed density rho_1 not initialized");
        if (this->grad.size() < 6) MSG_ERROR("Perturbed gradient not initialized");
        
        out_vec.push_back(this->rho[1]);
        out_vec.insert(out_vec.end(), this->grad.begin() + 3, this->grad.begin() + 6);
        
        // Include perturbed kinetic energy density for meta-GGA
        if (this->tau.size() >= 2) {
            out_vec.push_back(this->tau[1]);
        } else {
            MSG_ABORT("Perturbed kinetic energy density tau_1 not initialized for meta-GGA");
        }
    }
    return out_vec;
}

/** @brief Prepare input functions to xcfun
 *
 * Collects input densities and computes necessary gradients.
 *
 * Ordering of input:
 * inp_vec[0] = alpha_0
 * inp_vec[1] = beta_0
 * inp_vec[2] = alpha_1
 * inp_vec[3] = beta_1
 * ...
 */
// Why does this exist??? is never called in mrchem
void mGGA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != this->order) MSG_ERROR("Invalid input length");
    if (this->rho.size() > 0) MSG_ERROR("Density not empty");
    if (this->grad.size() > 0) MSG_ERROR("Gradient not empty");

    int n = 0;
    // For meta-GGA, expect alternating densities and tau: [rho_0, tau_0, rho_1, tau_1, ...]
    // For now, only support non-spin: [rho_0, tau_0] or [rho_0, tau_0, rho_1, tau_1]
    for (int i = 0; i < this->order; i++) {
        this->rho.push_back(inp_vec[n++]);
        if (n < static_cast<int>(inp_vec.size())) {
            this->tau.push_back(inp_vec[n++]);  // tau for each order
        }
    }

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

/** @brief Compute final output functions
 *
 * Combine the raw partial derivatives from xcfun into functional derivatives.
 *
 * For mGGA:
 * f_xc       : out[0] = inp[0]
 * df_xc/drho : out[1] = inp[1] - div(inp[2,3,4]) + inp[5]
 *
 * Note: inp[5] is v_tau (derivative w.r.t. kinetic energy density tau),
 * which is a local contribution (no divergence needed).
 */
mrcpp::FunctionTreeVector<3> mGGA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
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

    // Add v_tau contribution (local term, no divergence needed)
    if (inp_vec.size() > 5) {
        mrcpp::FunctionTree<3> &v_tau = mrcpp::get_func(inp_vec, 5);
        mrcpp::build_grid(*v_xc, v_tau);
        mrcpp::add(1.0, *v_xc, 1.0, v_tau);
    }

    // Collect output
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, v_xc));
    v_xc = nullptr;

    return out_vec;
}

} // namespace mrdft
