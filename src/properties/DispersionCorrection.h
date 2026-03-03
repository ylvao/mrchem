// /*
//  * MRChem, a numerical real-space code for molecular electronic structure
//  * calculations within the self-consistent field (SCF) approximations of quantum
//  * chemistry (Hartree-Fock and Density Functional Theory).
//  * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
//  *
//  * This file is part of MRChem.
//  *
//  * MRChem is free software: you can redistribute it and/or modify
//  * it under the terms of the GNU Lesser General Public License as published by
//  * the Free Software Foundation, either version 3 of the License, or
//  * (at your option) any later version.
//  *
//  * MRChem is distributed in the hope that it will be useful,
//  * but WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  * GNU Lesser General Public License for more details.
//  *
//  * You should have received a copy of the GNU Lesser General Public License
//  * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
//  *
//  * For information on the complete list of contributors to MRChem, see:
//  * <https://mrchem.readthedocs.io/>
//  */

#pragma once

#include <nlohmann/json.hpp>

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

class DispersionCorrection final {
public:
    explicit DispersionCorrection(int k = 1)
        : Nuclei nuc = this->getNuclei()),
          electronic(math_utils::init_nan(k, 3)) {}
        // // compute dispersion correction with external D3 library
        //     const Nuclei &nucs = this->getNuclei();
        //     int natom = static_cast<int>(nucs.size());
        //     std::vector<int> iz(natom);
        //     std::vector<double> coords(3 * natom);
        //     for (int i = 0; i < natom; ++i) {
        //         iz[i] = nucs[i].getElement().getZ();
        //         const auto &c = nucs[i].getCoord();
        //         coords[3 * i + 0] = c[0];
        //         coords[3 * i + 1] = c[1];
        //         coords[3 * i + 2] = c[2];
        //     }

        //     // parameters for the C wrapper defined in libd3_interface.h
        //     double energy_d3;
        //     std::vector<double> gradient(3 * natom);
        //     // version selects the parametrisation used by the D3 library.  a
        //     // value of 0 means "no version" and causes setfuncpar() to leave
        //     // all coefficients zero – hence the dispersion energy/gradients are
        //     // identically zero.  the Python examples use 4 (D3BJ variant) so we
        //     // follow suit here; the value could later be exposed as a user
        //     // option if desired.
        //     int version = 3;              // non‑zero -> use real parameters, 4 = BJ damping, zero damping = 3
        //     int tz = 0;                    // no tz scaling
        //     const char funcname[] = "pbe";  // currently hardcoded to PBE

         // call the Fortran-C wrapper
         // note: the interface uses value semantics for the scalar
         // integers (see updated libd3_interface.h), so pass them directly
         wrapper(natom,
                 coords.data(),      // array natoms-by-3 in row-major order
                 iz.data(),          // atomic numbers
                 funcname,            // functional name string (null terminated)
                 version,            // functional version
                 tz,                 // zero/one for tz flag
                 &energy_d3,         // output energy
                 gradient.data());   // output gradients (3 x natoms)

        //     // F.setDispersionCorrection(energy_d3);
        //     // std::cout << "[D3] dispersion energy = " << energy_d3
        //     //           << " au (functional='" << funcname << "', version="
        //     //           << version << ")\n";
        //     // std::cout << "Vector elements: ";
        //     //     for (const auto& element : gradient) {
        //     //         std::cout << element << " ";
        //     //     }
        //     // std::cout << std::endl;


        // // F.setDispersionCorrection(3.14);  // replaced by D3 call above
        // E_n = F.trace(Phi_n, nucs);
}
}