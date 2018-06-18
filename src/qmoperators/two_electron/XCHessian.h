#pragma once

#include "QMPotential.h"
#include "mrdft/XCFunctional.h"

/** @class XCHessian

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
