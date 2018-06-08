#pragma once

#include "qmfunctions.h"

/** @file sad.h
 *
 * @brief Module for generating initial guess as superposition of atomic densities
 *
 * The sad namespace provides a single function in the public interface,
 * used to generate an initial guess from hydrogen functions and a superposition
 * of atomic densities.
 */

namespace mrchem {
class Molecule;
class Nuclei;

namespace initial_guess {
namespace sad {

OrbitalVector setup(double prec, const Molecule &mol, bool restricted, int zeta);

} //namespace sad
} //namespace initial_guess
} //namespace mrchem
