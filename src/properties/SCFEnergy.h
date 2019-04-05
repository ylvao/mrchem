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

#include "mrchem.h"

/** @class SCFEnergy
 *
 * @brief Simple POD container to hold the different contributions to the SCF energy.
 *
 */

namespace mrchem {

// clang-format off
class SCFEnergy final {
public:
    SCFEnergy() = default;
    SCFEnergy(double nuc, double el, double orb,
              double kin, double en, double ee,
              double xc, double x, double nex,
              double ext) :
        E_nuc(nuc), E_el(el), E_orb(orb), E_kin(kin), E_en(en),
        E_ee(ee), E_x(x), E_xc(xc), E_nex(nex), E_ext(ext) {}

    double getTotalEnergy() const { return this->E_nuc + this->E_el; }
    double getNuclearEnergy() const { return this->E_nuc; }
    double getElectronicEnergy() const { return this->E_el; }

    double getOrbitalEnergy() const { return this->E_orb; }
    double getKineticEnergy() const { return this->E_kin; }
    double getElectronNuclearEnergy() const { return this->E_en; }
    double getElectronElectronEnergy() const { return this->E_ee; }
    double getExchangeCorrelationEnergy() const { return this->E_xc; }
    double getExchangeEnergy() const { return this->E_x; }

    friend std::ostream& operator<<(std::ostream &o, const SCFEnergy &en) {
        auto E_au = en.E_nuc + en.E_el;
        auto prec = mrcpp::Printer::getPrecision();

        std::stringstream o_kin, o_en, o_ee, o_x, o_xc, o_xe, o_xn;
        o_kin << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_kin;
        o_en << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_en;
        o_ee << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_ee;
        o_x << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_x;
        o_xc << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_xc;
        o_xe << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_ext;
        o_xn << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_nex;

        std::stringstream o_el, o_nuc;
        o_el << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_el;
        o_nuc << std::setw(27) << std::setprecision(2*prec) << std::fixed << en.E_nuc;

        std::stringstream o_au, o_ev, o_kj, o_kcal;
        o_au << std::setw(27) << std::setprecision(2*prec) << std::scientific << E_au;
        o_ev << std::setw(27) << std::setprecision(2*prec) << std::scientific << E_au * PHYSCONST::eV;
        o_kj << std::setw(27) << std::setprecision(2*prec) << std::scientific << E_au * PHYSCONST::kJ;
        o_kcal << std::setw(27) << std::setprecision(2*prec) << std::scientific << E_au * PHYSCONST::kcal;

        o << "============================================================" << std::endl;
        o << "                         SCF energy                         " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "        Kinetic :          (au) " << o_kin.str()              << std::endl;
        o << "            E-N :          (au) " << o_en.str()               << std::endl;
        o << "        Coulomb :          (au) " << o_ee.str()               << std::endl;
        o << "       Exchange :          (au) " << o_x.str()                << std::endl;
        o << "            X-C :          (au) " << o_xc.str()               << std::endl;
        o << " Ext. field (E) :          (au) " << o_xe.str()               << std::endl;
        o << " Ext. field (N) :          (au) " << o_xn.str()               << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "     Electronic :          (au) " << o_el.str()               << std::endl;
        o << "        Nuclear :          (au) " << o_nuc.str()              << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Total energy :          (au) "  << o_au.str()              << std::endl;
        o << "                :    (kcal/mol) "  << o_kcal.str()            << std::endl;
        o << "                :      (kJ/mol) "  << o_kj.str()              << std::endl;
        o << "                :          (eV) "  << o_ev.str()              << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;

        return o;
    }

private:
    double E_nuc{0.0};
    double E_el{0.0};

    double E_orb{0.0};
    double E_kin{0.0};
    double E_en{0.0};
    double E_ee{0.0};
    double E_x{0.0};
    double E_xc{0.0};
    double E_nex{0.0};
    double E_ext{0.0};
};
// clang-format on

} // namespace mrchem
