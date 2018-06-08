#pragma once

#include "qmfunctions.h"

/** @file core.h
 *
 * @brief Module for generating initial guess of hydrogen functions
 *
 * The initial_guess::core namespace provides a single function in the public
 * interface, used to set up an initial guess of hydrogen eigenfunctions.
 * The initial guess requires no external input, but is not very accurate.
 */

namespace mrchem {
class Molecule;
class Nuclei;

namespace initial_guess {
namespace core {

OrbitalVector setup(double prec, const Molecule &mol, bool restricted, int zeta);
OrbitalVector project_ao(double prec, const Nuclei &nucs, int spin, int zeta);

} //namespace core
} //namespace initial_guess
} //namespace mrchem
