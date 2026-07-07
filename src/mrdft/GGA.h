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

#pragma once

#include <XCFun/xcfun.h>

#include "Functional.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"

#include "xc_utils.h"


namespace mrdft {
class GGA final : public Functional {
public:
    GGA(int k, XClib_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d)
        : Functional(k, f)
        , derivative(std::move(d)) {
    xc_mask = xc_utils::build_output_mask(false, false, this->order);
    d_mask = xc_utils::build_density_mask(false, false, this->order);
    }
    ~GGA() override = default;

    bool isSpin() const override { return false; }
    bool isGGA() const override { return true; }
    bool isMetaGGA() const override { return false; }
    int numIn() const override { return 4; }
    int numOut() const override { return xclib->getnOut(); }
    int densityChannels() const override { return 1; }
    bool usesGradients() const override { return true; }

private:
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho;
    mrcpp::FunctionTreeVector<3> grad;

    int getCtrOutputLength() const override { return 5; }

    /** @brief Clear internal functions
     *
     * Ownership of densities is outside MRDFT -> clear
     * Ownership of gradients is inside MRDFT -> free
     */
    void clear() {
        mrcpp::clear(this->rho, false);
        mrcpp::clear(this->grad, true);
    }
};

} // namespace mrdft
