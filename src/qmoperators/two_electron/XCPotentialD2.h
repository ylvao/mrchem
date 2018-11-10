#pragma once

#include "qmoperators/one_electron/QMPotential.h"
#include "mrdft/XCFunctional.h"

/** @class XCPotential
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

class XCPotentialD2 final : public XCPotential {
public:
    XCPotentialD2(mrdft::XCFunctional *F,
                  OrbitalVector *Phi,
                  OrbitalVector *X,
                  OrbitalVector *Y);

protected:
    OrbitalVector *orbitals;                   ///< External set of orbitals used to build the density
    OrbitalVector *orbitals_x;                 ///< 1st external set of perturbed orbitals used to build the density
    OrbitalVector *orbitals_y;                 ///< 2nd external set of perturbed orbitals used to build the density
    mrcpp::FunctionTreeVector<3> potentials;   ///< XC Potential functions collected in a vector
    Density *pertDensity_t;                    ///< total first-order perturbed electronic density
    Density *pertDensity_a;                    ///< alpha first-order perturbed electronic density
    Density *pertDensity_b;                    ///< beta  first-order perturbed electronic density
  
    void setup(double prec);
    void clear();

    void setupDensity();
    void setupPotential(double prec);
    mrcpp::FunctionTree<3> &getPerturbedDensity(int spin);
    mrcpp::FunctionTree<3> &getPotential(int orbitalSpin, int densitySpin);

    Orbital apply(Orbital phi);

    //LUCA I wanted to include the following declarations in the cpp
    //file but i did not manage to get the syntax (copied from
    //density_utils.copp) right.
    int getPotentialIndex(int orbitalSpin, int densitySpin);
    void setupGroundStateDensity();
    void setupPerturbedDensity();

};

} //namespace mrchem
