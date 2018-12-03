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

OrbitalVector GroundStateHelmholtz::operator()(FockOperator &fock, const ComplexMatrix &F, OrbitalVector &Phi) const {
    Timer t_tot;
    Printer::printHeader(0, "Applying Helmholtz operators");
    int oldprec = Printer::setPrecision(5);

    Timer t_mat;
    ComplexVector lambda = F.diagonal();
    ComplexMatrix L = lambda.asDiagonal();
    OrbitalVector MPhi = orbital::rotate(L - F, Phi);
    t_mat.stop();
    Printer::printDouble(0, "Computing matrix part", t_mat.getWallTime(), 5);

    Printer::printSeparator(0, '-');
    println(0, " Orb    RealNorm   Nodes     ImagNorm   Nodes     Timing");
    Printer::printSeparator(0, '-');

    OrbitalVector out = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        Timer t_i;
        ComplexDouble mu_i = std::sqrt(-2.0 * lambda(i));
        if (std::abs(mu_i.imag()) > mrcpp::MachineZero) MSG_FATAL("Mu cannot be complex");

        Orbital psi_i = fock.potential()(Phi[i]);
        psi_i.add(1.0, MPhi[i]);
        psi_i.rescale(-1.0 / (2.0 * MATHCONST::pi));
        MPhi[i].free();

        mrcpp::HelmholtzOperator H_i(*MRA, mu_i.real(), this->build_prec);
        out[i] = apply(H_i, psi_i);
        psi_i.free();

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        t_i.stop();
        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(14) << rNorm);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(14) << iNorm);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, std::setw(14) << t_i.getWallTime() << std::endl);
    }

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
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
