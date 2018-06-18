#pragma once

#include "qmoperators/one_electron/QMPotential.h"
#include "qmfunctions/Density.h"

/** @class CoulombPotential
 *
 * @brief Coulomb potential defined by a particular electron density
 *
 * The Coulomb potential is computed by application of the Poisson operator
 * an electron density. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 */

namespace mrchem {

class CoulombPotential final : public QMPotential {
public:
    CoulombPotential(mrcpp::PoissonOperator *P, OrbitalVector *Phi = nullptr, OrbitalVector *Phi_1= nullptr, int order = 1);

    friend class CoulombOperator;

protected:
    Density density;                  ///< Ground-state electron density
    Density density_1;                ///< Perturbed electron density defining the potential
    OrbitalVector *orbitals;          ///< Unperturbed orbitals defining the ground-state electron density
    OrbitalVector *orbitals_1;        ///< Perturbed Orbitals defining the first order-electron density
    int order;                        ///< Order otf the potential (1=potential, 2=hessian, ....)
    mrcpp::PoissonOperator *poisson;  ///< Operator used to compute the potential

    Density &getDensity() { return this->density; }
    bool hasDensity() const { return (this->density.getSquareNorm() < 0.0) ? false : true; }
    bool hasDensity_1() const { return (this->density_1.getSquareNorm() < 0.0) ? false : true; }

    void setup(double prec);
    void clear();

    void setupDensity(double prec);
    void setupPotential(double prec);
    void setupHessian(double prec);
};

} //namespace mrchem
