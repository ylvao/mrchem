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

#include "SpinGGA.h"
#include "xc_utils.h"

namespace mrdft {

SpinGGA::SpinGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d)
        : Functional(k, f)
        , derivative(std::move(d)) {
    xc_mask = xc_utils::build_output_mask(false, true, this->order);
    d_mask = xc_utils::build_density_mask(false, true, this->order);
}

/** @brief Clear internal functions
 *
 * Ownership of densities is outside MRDFT -> clear
 * Ownership of gradients is inside MRDFT -> free
 */
void SpinGGA::clear() {
    mrcpp::clear(this->rho_a, false);
    mrcpp::clear(this->rho_b, false);
    mrcpp::clear(this->grad_a, true);
    mrcpp::clear(this->grad_b, true);
}

/** @brief Number of function involved in contraction step */
int SpinGGA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 8;
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Collect input functions to xcfun evaluation step
 *
 * For SpinGGA : [alpha_0, beta_0, grad(alpha_0), grad(beta_0)]
 */
mrcpp::FunctionTreeVector<3> SpinGGA::setupXCInput() {
    if (this->rho_a.size() < 1) MSG_ERROR("Alpha density not initialized");
    if (this->rho_b.size() < 1) MSG_ERROR("Beta density not initialized");
    if (this->grad_a.size() < 3) MSG_ERROR("Alpha gradient not initialized");
    if (this->grad_b.size() < 3) MSG_ERROR("Beta gradient not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho_a[0]);
    out_vec.push_back(this->rho_b[0]);
    out_vec.insert(out_vec.end(), this->grad_a.begin(), this->grad_a.begin() + 3);
    out_vec.insert(out_vec.end(), this->grad_b.begin(), this->grad_b.begin() + 3);
    return out_vec;
}

/** @brief Collect input functions to contraction step
 *
 * For SpinGGA:
 * Ground State: No contraction, empty vector
 * Linear Response: [alpha_1, beta_1, grad(alpha_1), grad(beta_1)]
 * Higher Response: NOT_IMPLEMENTED
 */
mrcpp::FunctionTreeVector<3> SpinGGA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (this->order == 2) {
        out_vec.push_back(this->rho_a[1]);
        out_vec.push_back(this->rho_b[1]);
        out_vec.insert(out_vec.end(), this->grad_a.begin() + 3, this->grad_a.begin() + 6);
        out_vec.insert(out_vec.end(), this->grad_b.begin() + 3, this->grad_b.begin() + 6);
    }
    return out_vec;
}

} // namespace mrdft
