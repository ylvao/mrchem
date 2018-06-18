#pragma once

#include "QMPotential.h"
#include "mrdft/XCFunctional.h"

/** @class XCHessian
 *
 * @brief Exchange-Correlation potential defined by a particular (spin) density
 *
 * The XC potential is computed by mapping of the density through a XC functional,
 * provided by the XCFun library. There are two ways of defining the density:
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
 *
 * LDA and GGA functionals are supported as well as two different ways to compute
 * the XC potentials: either with explicit derivatives or gamma-type derivatives.
 */

namespace mrchem {

class XCHessian final : public QMPotential {
public:
    XCHessian(mrdft::XCFunctional *F, OrbitalVector *Phi = nullptr);

    friend class XCOperator;

protected:
    OrbitalVector *orbitals;                   ///< External set of orbitals used to build the density
    mrdft::XCFunctional *functional;           ///< External XC functional to be used

    mrcpp::FunctionTreeVector<3> hessians ;   ///< XC Hessian components collected in a vector

    Density &getDensity(int spin);
    mrcpp::FunctionTree<3> &getHessianComponent(int spin_den, int spin_orb);

    void setup(double prec);
    void clear();

    void setupDensity();
    void setupHessian(double prec);

    Orbital apply(Orbital phi);
};

} //namespace mrchem
