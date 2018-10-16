/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "mrchem.h"

namespace mrchem {

class Polarizability final {
public:
    Polarizability(double w = 0.0, const double *o = 0, bool v = false) {
        this->frequency = w;
        this->origin[0] = (o != 0) ? o[0] : 0.0;
        this->origin[1] = (o != 0) ? o[1] : 0.0;
        this->origin[2] = (o != 0) ? o[2] : 0.0;
        this->velocity = v;
        this->tensor = DoubleMatrix::Zero(3,3);
    }
    ~Polarizability() { }

    double getFrequency() const { return this->frequency; }
    bool getVelocityGauge() const { return this->velocity; }
    const double* getGaugeOrigin() const { return this->origin; }
    DoubleMatrix &get() { return this->tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Polarizability &pol) {
        DoubleMatrix tens = pol.tensor;

        double w_au = 0.0;  // Only static Polarizability
        double isoPolarau = tens.trace()/3.0;

        double isoPolarsi = isoPolarau * 0.0; // Luca: FIX THIS

        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Polarizability tensor                    "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Frequency:       (au)       " << std::setw(30) << w_au      <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Isotropic average(au)       " << std::setw(30) << isoPolarau<<std::endl;
        o<<" Isotropic average(SI)       TO BE FIXED                    "<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Tensor --------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<tens(0,0)<<std::setw(20)<<tens(1,0)<<std::setw(20)<<tens(2,0)<<std::endl;
        o<<std::setw(19)<<tens(1,0)<<std::setw(20)<<tens(1,1)<<std::setw(20)<<tens(2,1)<<std::endl;
        o<<std::setw(19)<<tens(2,0)<<std::setw(20)<<tens(1,2)<<std::setw(20)<<tens(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    double frequency;
    double origin[3];
    bool velocity;
    DoubleMatrix tensor;
};

} //namespace mrchem
