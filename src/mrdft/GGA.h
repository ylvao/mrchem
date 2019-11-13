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

#pragma once

#include "Functional.h"

/**
 *  @class XCFunctional
 *  @brief Compute XC functional with XCFun
 *
 *  Interface class for the XCFun library
 *
 * Depending on the mode chosen, xcfun needs either the gamma
 * functions or the explicit gradients. The first mode is possibly
 * more efficient (fewer functions to compute/handle), whereas the
 * other is simpler to implement. We keep both options open and
 * compute the gradient invariants if and when necessary.
 *
 * The XCFunctional keeps track of the density grid, whose size is
 * initially defined through the interface function buildGrid().
 * The grid is kept _fixed_ for all internal calculations
 * within the module.
 *
 * Typical usage within one SCF cycle:
 *
 * 1) getDensity() and compute density on given grid
 * 2) setup()
 * 3) evaluate()
 * 4) calcEnergy()
 * 5) calcPotential()
 * 7) clear()
 *
 */

namespace mrdft {

class GGA final : public Functional {
public:
    GGA(int k, std::unique_ptr<xc_functional> &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d);
    ~GGA() override = default;

    bool isSpin() const override { return false; }

private:
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho;
    mrcpp::FunctionTreeVector<3> grad;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 5; }

    void clear() override;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() override;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
