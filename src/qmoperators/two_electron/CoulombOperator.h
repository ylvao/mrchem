#pragma once

#include "qmoperators/RankZeroTensorOperator.h"
#include "CoulombPotential.h"
#include "CoulombHessian.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class CoulombOperator final : public RankZeroTensorOperator {
public:
    CoulombOperator(mrcpp::PoissonOperator *P,
                    OrbitalVector *Phi = nullptr) {
        this->potential = new CoulombPotential(P, Phi);
        this->hessian = nullptr;
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }

    ~CoulombOperator() {
        if (this->potential != nullptr) delete this->potential;
        if (this->hessian != nullptr) delete this->hessian;
    }
    
    CoulombOperator(mrcpp::PoissonOperator *P,
                    OrbitalVector *Phi,
                    OrbitalVector *Phi_x) {
        this->potential = nullptr;
        this->hessian = new CoulombHessian(P, Phi, Phi_x);
        RankZeroTensorOperator &J = (*this);
        J = *this->hessian;
    }
    
    CoulombOperator(mrcpp::PoissonOperator *P,
                    OrbitalVector *Phi,
                    OrbitalVector *Phi_x,
                    OrbitalVector *Phi_y) {
        this->potential = nullptr;
        this->hessian = new CoulombHessian(P, Phi, Phi_x, Phi_y);
        RankZeroTensorOperator &J = (*this);
        J = *this->hessian;
    }

    Density &getDensity() {
        if (potential != nullptr) return this->potential->getDensity();
        if (hessian != nullptr) return this->hessian->getDensity();
        MSG_FATAL("Coulomb operator not properly initialized");
    }
        
    ComplexDouble trace(OrbitalVector &Phi) { return 0.5*RankZeroTensorOperator::trace(Phi); }

 protected:
    CoulombPotential *potential;
    CoulombHessian *hessian;
};

} //namespace mrchem
