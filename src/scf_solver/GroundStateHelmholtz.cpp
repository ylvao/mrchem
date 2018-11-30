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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "GroundStateHelmholtz.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

OrbitalVector GroundStateHelmholtz::operator()(FockOperator &fock, ComplexMatrix &F, OrbitalVector &inp) {
    Printer::printHeader(0, "Applying Helmholtz operators");
    println(0, " Orb    RealNorm   Nodes     ImagNorm   Nodes     Timing");
    Printer::printSeparator(0, '-');
    int oldprec = Printer::setPrecision(5);

    OrbitalVector out = orbital::param_copy(inp);

    QMFunctionVector phi_vec;
    for (int i = 0; i < inp.size(); i++) phi_vec.push_back(inp[i]);

    Timer tottimer;
    for (int i = 0; i < inp.size(); i++) {
        Timer timer;
        ComplexDouble mu = std::sqrt(-2.0 * F(i, i));
        if (std::abs(mu.imag()) > mrcpp::MachineZero) MSG_FATAL("Mu cannot be complex");

        ComplexVector c = -1.0 * F.row(i);
        c(i) = 0.0;

        Orbital chi = inp[i].paramCopy();
        qmfunction::linear_combination(chi, c, phi_vec, this->apply_prec);

        Orbital psi = fock.potential()(inp[i]);
        psi.add(1.0, chi);
        psi.rescale(-1.0 / (2.0 * MATHCONST::pi));
        chi.free();

        mrcpp::HelmholtzOperator H(*MRA, mu.real(), this->build_prec);
        out[i] = apply(H, psi);
        psi.free();

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        timer.stop();
        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(14) << rNorm);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(14) << iNorm);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, std::setw(14) << timer.getWallTime() << std::endl);
    }

    tottimer.stop();
    Printer::printFooter(0, tottimer, 2);
    Printer::setPrecision(oldprec);
    return out;
}

Orbital GroundStateHelmholtz::apply(mrcpp::HelmholtzOperator &H, Orbital &inp) const {
    Orbital out = inp.paramCopy();
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(this->apply_prec, out.real(), H, inp.real());
    }
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(this->apply_prec, out.imag(), H, inp.imag());
        if (inp.conjugate()) out.imag().rescale(-1.0);
    }
    return out;
}
} // namespace mrchem
