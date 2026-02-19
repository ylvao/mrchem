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

// #pragma once

// #include <nlohmann/json.hpp>

// #include "mrchem.h"

// #include "utils/math_utils.h"
// #include "utils/print_utils.h"

// namespace mrchem {

// class DispersionCorrection final {
// public:

// /**
//   * @brief Helper function to find a C6 coefficient for a pair of atoms
//   * based on their atomic numbers and coordination numbers.
//   */
//  void Molecule::printDispersionForMolecule() {
//      const auto& nuclei = mol.getNuclei();
//      int nAtoms = mol.getNNuclei();

//      std::cout << "--- Dispersion Coefficient Lookup for Molecule ---" << std::endl;

//      // Iterate over all unique pairs of atoms (i, j)
//      for (int i = 0; i < nAtoms; ++i) {
//          for (int j = i + 1; j < nAtoms; ++j) {
//              int z1 = nuclei[i].getElement().getZ(); // Get Z from Molecule
//              int z2 = nuclei[j].getElement().getZ(); // Get Z from Molecule

//              bool pairFound = false;

//              // Search the parsed DISPERSION_DATA array
//              for (int k = 0; k < DISPERSION_DATA_SIZE; ++k) {
//                  const auto& entry = DISPERSION_DATA[k];

//                  // Check for matching atomic numbers (order independent)
//                  if ((entry.z1 == z1 && entry.z2 == z2) || (entry.z1 == z2 && entry.z2 == z1)) {
//                      std::cout << "Pair " << i << "(" << nuclei[i].getElement().getSymbol() << ") - "
//                                << j << "(" << nuclei[j].getElement().getSymbol() << "): "
//                                << "C6 = " << entry.c6 << " (at ref CN1=" << entry.cn1 
//                                << ", CN2=" << entry.cn2 << ")" << std::endl;
//                      pairFound = true;
//                  }
//              }

//              if (!pairFound) {
//                  std::cout << "No D3 parameters found for pair Z=(" << z1 << "," << z2 << ")" << std::endl;
//              }
//          }
//      }
//  }


// }
// }