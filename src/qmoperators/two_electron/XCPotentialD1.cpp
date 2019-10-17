#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "XCPotentialD1.h"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

XCPotentialD1::XCPotentialD1(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi, bool mpi_shared)
        : XCPotential(F, Phi, mpi_shared) {
    densities.push_back(Density(false)); // rho_0 total
    densities.push_back(Density(false)); // rho_0 alpha
    densities.push_back(Density(false)); // rho_0 beta
}

mrcpp::FunctionTreeVector<3> XCPotentialD1::setupDensities(double prec, mrcpp::FunctionTree<3> &grid) {
    mrcpp::FunctionTreeVector<3> dens_vec;
    if (not this->mrdft->functional().isSpin()) {
        { // Unperturbed total density
            Timer timer;
            Density &rho = getDensity(DENSITY::DensityType::Total, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DENSITY::DensityType::Total);
            }
            print_utils::qmfunction(2, "XC rho", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
    } else {
        { // Unperturbed alpha density
            Timer timer;
            Density &rho = getDensity(DENSITY::DensityType::Alpha, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DENSITY::DensityType::Alpha);
            }
            print_utils::qmfunction(2, "XC rho (alpha)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Unperturbed beta density
            Timer timer;
            Density &rho = getDensity(DENSITY::DensityType::Beta, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DENSITY::DensityType::Beta);
            }
            print_utils::qmfunction(2, "XC rho (beta)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
    }
    return dens_vec;
}

} // namespace mrchem
