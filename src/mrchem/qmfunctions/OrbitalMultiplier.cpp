#include "OrbitalMultiplier.h"
#include "Orbital.h"
#include "Potential.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

OrbitalMultiplier::OrbitalMultiplier(double prec)
    : add(prec, MRA->getMaxScale()),
      mult(prec, MRA->getMaxScale()),
      grid(MRA->getMaxScale()) {
}

void OrbitalMultiplier::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void OrbitalMultiplier::operator()(Orbital &Vphi, Potential &V, Orbital &phi) {
    calcRealPart(Vphi, V, phi);
    calcImagPart(Vphi, V, phi, false);
}

void OrbitalMultiplier::adjoint(Orbital &Vphi, Potential &V, Orbital &phi) {
    calcRealPart(Vphi, V, phi);
    calcImagPart(Vphi, V, phi, true);
}

// phi_ab = c * phi_a * phi_b
void OrbitalMultiplier::operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    calcRealPart(phi_ab, c, phi_a, phi_b);
    calcImagPart(phi_ab, c, phi_a, phi_b, false);
}

// phi_ab = c * phi_a^dag * phi_b
void OrbitalMultiplier::adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    calcRealPart(phi_ab, c, phi_a, phi_b);
    calcImagPart(phi_ab, c, phi_a, phi_b, true);
}

void OrbitalMultiplier::calcRealPart(Orbital &Vphi,
                                     Potential &V,
                                     Orbital &phi) {
    if (Vphi.hasReal()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> vec;
    if (V.hasReal() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.re());
        this->mult(*tree, 1.0, V.re(), phi.re(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.im());
        this->mult(*tree, 1.0, V.im(), phi.im(), 0);
        vec.push_back(-1.0, tree);
    }
    if (vec.size() == 1) {
        Vphi.setReal(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        Vphi.allocReal();
        this->grid(Vphi.re(), vec);
        this->add(Vphi.re(), vec, 0);
        vec.clear(true);
    }
}

void OrbitalMultiplier::calcImagPart(Orbital &Vphi,
                                     Potential &V,
                                     Orbital &phi,
                                     bool adjoint) {
    if (Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> vec;
    if (V.hasReal() and phi.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.im());
        this->mult(*tree, 1.0, V.re(), phi.im(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.re());
        this->mult(*tree, 1.0, V.im(), phi.re(), 0);
        if (adjoint) {
            vec.push_back(-1.0, tree);
        } else {
            vec.push_back(1.0, tree);
        }
    }
    if (vec.size() == 1) {
        Vphi.setImag(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        Vphi.allocImag();
        this->grid(Vphi.im(), vec);
        this->add(Vphi.im(), vec, 0);
        vec.clear(true);
    }
}

void OrbitalMultiplier::calcRealPart(Orbital &phi_ab,
                                     double c,
                                     Orbital &phi_a,
                                     Orbital &phi_b) {
    if (phi_ab.hasReal()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.re());
        this->grid(*tree, phi_b.re());
        this->mult(*tree, c, phi_a.re(), phi_b.re(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.im());
        this->grid(*tree, phi_b.im());
        this->mult(*tree, c, phi_a.im(), phi_b.im(), 0);
        vec.push_back(-1.0, tree);
    }
    if (vec.size() == 1) {
        phi_ab.setReal(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        phi_ab.allocReal();
        this->grid(phi_ab.re(), vec);
        this->add(phi_ab.re(), vec, 0);
        vec.clear(true);
    }
}

void OrbitalMultiplier::calcImagPart(Orbital &phi_ab,
                                     double c,
                                     Orbital &phi_a,
                                     Orbital &phi_b,
                                     bool adjoint) {
    if (phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.re());
        this->grid(*tree, phi_b.im());
        this->mult(*tree, c, phi_a.re(), phi_b.im(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.im());
        this->grid(*tree, phi_b.re());
        this->mult(*tree, c, phi_a.im(), phi_b.re(), 0);
        if (adjoint) {
            vec.push_back(-1.0, tree);
        } else {
            vec.push_back(1.0, tree);
        }
    }
    if (vec.size() == 1) {
        phi_ab.setImag(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        phi_ab.allocImag();
        this->grid(phi_ab.im(), vec);
        this->add(phi_ab.im(), vec, 0);
        vec.clear(true);
    }
}
