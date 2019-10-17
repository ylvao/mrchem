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

#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <XCFun/xcfun.h>

namespace mrdft {

class Functional {
public:
    Functional(int k, std::unique_ptr<xc_functional> &f)
            : order(k)
            , xcfun(std::move(f)) {}
    virtual ~Functional() = default;

    virtual bool isSpin() const = 0;
    bool isLDA() const { return (not(isGGA() or isMetaGGA())); }
    bool isGGA() const { return xc_is_gga(*xcfun); }
    bool isMetaGGA() const { return xc_is_metagga(*xcfun); }
    bool isHybrid() const { return (std::abs(amountEXX()) > 1.0e-10); }
    double amountEXX() const {
        double exx = 0.0;
        xc_get(*xcfun, "exx", &exx);
        return exx;
    }

    friend class MRDFT;

protected:
    const int order;
    Eigen::VectorXi d_mask;
    Eigen::MatrixXi xc_mask;
    std::unique_ptr<xc_functional> xcfun;

    int getXCInputLength() const { return xc_input_length(*xcfun); }
    int getXCOutputLength() const { return xc_output_length(*xcfun); }
    virtual int getCtrInputLength() const = 0;
    virtual int getCtrOutputLength() const = 0;

    Eigen::MatrixXd evaluate(Eigen::MatrixXd &inp) const;
    Eigen::MatrixXd contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;

    virtual void clear() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() = 0;

    virtual void preprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
    virtual mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
};

} // namespace mrdft
