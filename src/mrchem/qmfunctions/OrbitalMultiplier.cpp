#include "OrbitalMultiplier.h"
#include "Orbital.h"
#include "Potential.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

OrbitalMultiplier::OrbitalMultiplier(double prec, int max_scale)
    : add(prec, max_scale),
      mult(prec, max_scale),
      grid(max_scale) {
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

void OrbitalMultiplier::calcRealPart(Orbital &Vphi, Potential &V, Orbital &phi) {
    if (Vphi.hasReal()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> vec;
    if (V.hasReal() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.real());
        this->mult(*tree, 1.0, V.real(), phi.real(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.imag());
        this->mult(*tree, 1.0, V.imag(), phi.imag(), 0);
        vec.push_back(-1.0, tree);
    }
    if (vec.size() == 1) {
        Vphi.setReal(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        Vphi.allocReal();
        this->grid(Vphi.real(), vec);
        this->add(Vphi.real(), vec, 0);
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
        this->grid(*tree, phi.imag());
        this->mult(*tree, 1.0, V.real(), phi.imag(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi.real());
        this->mult(*tree, 1.0, V.imag(), phi.real(), 0);
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
        this->grid(Vphi.imag(), vec);
        this->add(Vphi.imag(), vec, 0);
        vec.clear(true);
    }
}

void OrbitalMultiplier::calcRealPart(Orbital &phi_ab,
                                     double c,
                                     Orbital &phi_a,
                                     Orbital &phi_b) {
    if (phi_ab.hasReal()) MSG_ERROR("Orbital not empty");

    // sanity check spin
    if (phi_ab.getSpin() != phi_a.getSpin()) MSG_FATAL("Mixing spins");
    if (phi_ab.getSpin() != phi_b.getSpin()) MSG_FATAL("Mixing spins");

    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.real());
        this->grid(*tree, phi_b.real());
        this->mult(*tree, c, phi_a.real(), phi_b.real(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.imag());
        this->grid(*tree, phi_b.imag());
        this->mult(*tree, c, phi_a.imag(), phi_b.imag(), 0);
        vec.push_back(-1.0, tree);
    }
    if (vec.size() == 1) {
        phi_ab.setReal(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        phi_ab.allocReal();
        this->grid(phi_ab.real(), vec);
        this->add(phi_ab.real(), vec, 0);
        vec.clear(true);
    }
}

void OrbitalMultiplier::calcImagPart(Orbital &phi_ab,
                                     double c,
                                     Orbital &phi_a,
                                     Orbital &phi_b,
                                     bool adjoint) {
    if (phi_ab.hasImag()) MSG_ERROR("Orbital not empty");

    // sanity check spin
    if (phi_ab.getSpin() != phi_a.getSpin()) MSG_FATAL("Mixing spins");
    if (phi_ab.getSpin() != phi_b.getSpin()) MSG_FATAL("Mixing spins");

    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.real());
        this->grid(*tree, phi_b.imag());
        this->mult(*tree, c, phi_a.real(), phi_b.imag(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        this->grid(*tree, phi_a.imag());
        this->grid(*tree, phi_b.real());
        this->mult(*tree, c, phi_a.imag(), phi_b.real(), 0);
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
        this->grid(phi_ab.imag(), vec);
        this->add(phi_ab.imag(), vec, 0);
        vec.clear(true);
    }
}
