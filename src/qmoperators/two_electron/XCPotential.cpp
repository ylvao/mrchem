#include "XCPotential.h"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Compute electron density
 *
 * The density is computed on the grid provided by the MRDFT module. The grid
 * is kept as is, e.i. no additional refinement at this point, since the grid
 * size is determined inside the module.
 */
void XCPotential::setupDensity(double prec) {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        FunctionTree<3> &func_a = this->getDensity(DENSITY::Alpha);
        Density rho_a(false);
        rho_a.setReal(&func_a);
        density::compute(prec, rho_a, Phi, DENSITY::Alpha);
        print_utils::qmfunction(1, "XC alpha density", rho_a, time_a);
        rho_a.setReal(nullptr);

        Timer time_b;
        FunctionTree<3> &func_b = this->getDensity(DENSITY::Beta);
        Density rho_b(false);
        rho_b.setReal(&func_b);
        density::compute(prec, rho_b, Phi, DENSITY::Beta);
        print_utils::qmfunction(1, "XC beta density", rho_b, time_b);
        rho_b.setReal(nullptr);

        // Extend to union grid
        while (mrcpp::refine_grid(func_a, func_b)) {}
        while (mrcpp::refine_grid(func_b, func_a)) {}
    } else {
        Timer time_t;
        FunctionTree<3> &func_t = this->getDensity(DENSITY::Total);
        Density rho_t(false);
        rho_t.setReal(&func_t);
        density::compute(prec, rho_t, Phi, DENSITY::Total);
        print_utils::qmfunction(1, "XC total density", rho_t, time_t);
        rho_t.setReal(nullptr);
    }
}

mrcpp::FunctionTree<3> &XCPotential::getDensity(int spin) {
    if (spin == DENSITY::Total) return this->functional->getDensity(mrdft::DensityType::Total);
    if (spin == DENSITY::Alpha) return this->functional->getDensity(mrdft::DensityType::Alpha);
    if (spin == DENSITY::Beta) return this->functional->getDensity(mrdft::DensityType::Beta);
    MSG_ABORT("Invalid density type");
}

// NOTE AFTER DISCUSSION WITH STIG: Need to move stuff that is
// iteration-independent out of the response loop, so that all required
// functions, which only depend on the GS density are computed
// once. Comment: still the grid for rho_1 is borrowed from rho_0 and
// rho_0 should still be available.
} // namespace mrchem
