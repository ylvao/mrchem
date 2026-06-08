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
#include "Factory.h" // only to call Factory::libxc

namespace mrdft {

class mGGA final : public Functional {
public:
    mGGA(int k, XC_p &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> &d);
    ~mGGA() override = default;

    bool isSpin() const override { return false; }
    bool isGGA() const override { return false; }
    bool isMetaGGA() const override { return true; }
    int numIn() const override { return 5; }
    int numOut() const override {
        if (Factory::libxc) {
            // LibXC mGGA outputs: exc, vxc, vsigma, v_laplacian (optional), v_tau.
            // We expose:
            //  row 0: f_xc * rho
            //  row 1: v_rho
            //  row 2-4: gradient contributions (x, y, z)
            //  row 5: v_tau
            return 6;
        } else {
            return xcfun_output_length(xcfun.get());
        }
    }
private:
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho;
    mrcpp::FunctionTreeVector<3> grad;
    mrcpp::FunctionTreeVector<3> tau;  ///< Kinetic energy density for meta-GGA

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 5; }

    void clear() override;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() override;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
