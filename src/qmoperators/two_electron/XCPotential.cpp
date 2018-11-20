#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "XCPotential.h"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

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
        FunctionTree<3> &func_a = this->functional->getDensity(mrdft::DensityType::Alpha);
        Density rho_a(false);
        rho_a.set(NUMBER::Real, &func_a);
        density::compute(prec, rho_a, Phi, DENSITY::Alpha);
        rho_a.set(NUMBER::Real, nullptr);
        time_a.stop();
        Printer::printTree(0, "XC alpha density", func_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        FunctionTree<3> &func_b = this->functional->getDensity(mrdft::DensityType::Beta);
        Density rho_b(false);
        rho_b.set(NUMBER::Real, &func_b);
        density::compute(prec, rho_b, Phi, DENSITY::Beta);
        rho_b.set(NUMBER::Real, nullptr);
        time_b.stop();
        Printer::printTree(0, "XC beta density", func_b.getNNodes(), time_b.getWallTime());

        // Extend to union grid
        while (mrcpp::refine_grid(func_a, func_b)) {}
        while (mrcpp::refine_grid(func_b, func_a)) {}
    } else {
        Timer time_t;
        FunctionTree<3> &func_t = this->functional->getDensity(mrdft::DensityType::Total);
        Density rho_t(false);
        rho_t.set(NUMBER::Real, &func_t);
        density::compute(prec, rho_t, Phi, DENSITY::Total);
        rho_t.set(NUMBER::Real, nullptr);
        time_t.stop();
        Printer::printTree(0, "XC total density", func_t.getNNodes(), time_t.getWallTime());
    }
}

}
