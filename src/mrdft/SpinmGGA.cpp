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

#include "SpinmGGA.h"
#include "xc_utils.h"

namespace mrdft {

SpinmGGA::SpinmGGA(int k, XC_p &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> &d)
        : Functional(k, f)
        , derivative(d) {
    xc_mask = xc_utils::build_output_mask(false, true, this->order);
    d_mask = xc_utils::build_density_mask(false, true, this->order);
}

void SpinmGGA::clear() {
    mrcpp::clear(this->rho_a, false);
    mrcpp::clear(this->rho_b, false);
    mrcpp::clear(this->grad_a, true);
    mrcpp::clear(this->grad_b, true);
    mrcpp::clear(this->tau_a, true);
    mrcpp::clear(this->tau_b, true);
}

int SpinmGGA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 10;  // rho_a_1, rho_b_1, grad_a_1×3, grad_b_1×3, tau_a_1, tau_b_1
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

mrcpp::FunctionTreeVector<3> SpinmGGA::setupXCInput() {
    if (this->rho_a.size() < 1) MSG_ERROR("Alpha density not initialized");
    if (this->rho_b.size() < 1) MSG_ERROR("Beta density not initialized");
    if (this->grad_a.size() < 3) MSG_ERROR("Alpha gradient not initialized");
    if (this->grad_b.size() < 3) MSG_ERROR("Beta gradient not initialized");
    if (this->tau_a.size() < 1) MSG_ERROR("Alpha kinetic energy density not initialized");
    if (this->tau_b.size() < 1) MSG_ERROR("Beta kinetic energy density not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho_a[0]);
    out_vec.push_back(this->rho_b[0]);
    out_vec.insert(out_vec.end(), this->grad_a.begin(), this->grad_a.begin() + 3);
    out_vec.insert(out_vec.end(), this->grad_b.begin(), this->grad_b.begin() + 3);
    out_vec.push_back(this->tau_a[0]);
    out_vec.push_back(this->tau_b[0]);
    return out_vec;
}

mrcpp::FunctionTreeVector<3> SpinmGGA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (this->order == 2) {
        if (this->rho_a.size() < 2) MSG_ERROR("Perturbed alpha density not initialized");
        if (this->rho_b.size() < 2) MSG_ERROR("Perturbed beta density not initialized");
        
        out_vec.push_back(this->rho_a[1]);
        out_vec.push_back(this->rho_b[1]);
        out_vec.insert(out_vec.end(), this->grad_a.begin() + 3, this->grad_a.begin() + 6);
        out_vec.insert(out_vec.end(), this->grad_b.begin() + 3, this->grad_b.begin() + 6);
        
        // Include perturbed kinetic energy densities
        if (this->tau_a.size() >= 2 && this->tau_b.size() >= 2) {
            out_vec.push_back(this->tau_a[1]);
            out_vec.push_back(this->tau_b[1]);
        } else {
            MSG_ABORT("Perturbed kinetic energy densities not initialized for spin meta-GGA");
        }
    }
    return out_vec;
}

void SpinmGGA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != this->order * 2) MSG_ERROR("Invalid input length for spin meta-GGA");
    if (this->rho_a.size() > 0) MSG_ERROR("Alpha density not empty");
    if (this->rho_b.size() > 0) MSG_ERROR("Beta density not empty");
    if (this->tau_a.size() > 0) MSG_ERROR("Alpha tau not empty");
    if (this->tau_b.size() > 0) MSG_ERROR("Beta tau not empty");

    int n = 0;
    for (int i = 0; i < this->order; i++) {
        this->rho_a.push_back(inp_vec[n++]);
        this->rho_b.push_back(inp_vec[n++]);
        if (n < static_cast<int>(inp_vec.size())) {
            this->tau_a.push_back(inp_vec[n++]);
            this->tau_b.push_back(inp_vec[n++]);
        }
    }

    for (int i = 0; i < this->order; i++) {
        mrcpp::FunctionTreeVector<3> tmp_a = mrcpp::gradient(*this->derivative, mrcpp::get_func(this->rho_a, i));
        mrcpp::FunctionTreeVector<3> tmp_b = mrcpp::gradient(*this->derivative, mrcpp::get_func(this->rho_b, i));
        this->grad_a.insert(this->grad_a.end(), tmp_a.begin(), tmp_a.end());
        this->grad_b.insert(this->grad_b.end(), tmp_b.begin(), tmp_b.end());
    }
}

mrcpp::FunctionTreeVector<3> SpinmGGA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Energy density
    mrcpp::FunctionTree<3> &f_xc = mrcpp::get_func(inp_vec, 0);
    inp_vec[0] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Potential parts
    mrcpp::FunctionTree<3> &df_dr_a = mrcpp::get_func(inp_vec, 1);
    mrcpp::FunctionTree<3> &df_dr_b = mrcpp::get_func(inp_vec, 2);
    mrcpp::FunctionTreeVector<3> df_dg_a(inp_vec.begin() + 3, inp_vec.begin() + 6);
    mrcpp::FunctionTreeVector<3> df_dg_b(inp_vec.begin() + 6, inp_vec.begin() + 9);

    auto *tmp_a = new mrcpp::FunctionTree<3>(df_dr_a.getMRA());
    auto *tmp_b = new mrcpp::FunctionTree<3>(df_dr_b.getMRA());
    mrcpp::divergence(*tmp_a, *this->derivative, df_dg_a);
    mrcpp::divergence(*tmp_b, *this->derivative, df_dg_b);

    auto *v_xc_a = new mrcpp::FunctionTree<3>(df_dr_a.getMRA());
    auto *v_xc_b = new mrcpp::FunctionTree<3>(df_dr_b.getMRA());
    mrcpp::build_grid(*v_xc_a, df_dr_a);
    mrcpp::build_grid(*v_xc_a, *tmp_a);
    mrcpp::build_grid(*v_xc_b, df_dr_b);
    mrcpp::build_grid(*v_xc_b, *tmp_b);
    mrcpp::add(-1.0, *v_xc_a, 1.0, df_dr_a, -1.0, *tmp_a);
    mrcpp::add(-1.0, *v_xc_b, 1.0, df_dr_b, -1.0, *tmp_b);
    delete tmp_a;
    delete tmp_b;

    // Add v_tau contributions (local terms, no divergence needed)
    if (inp_vec.size() > 9) {
        mrcpp::FunctionTree<3> &v_tau_a = mrcpp::get_func(inp_vec, 9);
        mrcpp::FunctionTree<3> &v_tau_b = mrcpp::get_func(inp_vec, 10);
        mrcpp::build_grid(*v_xc_a, v_tau_a);
        mrcpp::build_grid(*v_xc_b, v_tau_b);
        mrcpp::add(1.0, *v_xc_a, 1.0, v_tau_a);
        mrcpp::add(1.0, *v_xc_b, 1.0, v_tau_b);
    }

    // Collect output
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, v_xc_a));
    out_vec.push_back(std::make_tuple(1.0, v_xc_b));

    return out_vec;
}

} // namespace mrdft
