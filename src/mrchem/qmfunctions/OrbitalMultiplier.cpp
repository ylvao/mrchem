#include "OrbitalMultiplier.h"
#include "Orbital.h"
#include "Potential.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

OrbitalMultiplier::OrbitalMultiplier(double pr)
    : add(*MRA, pr),
      mult(*MRA, pr),
      grid(*MRA) {
}

void OrbitalMultiplier::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void OrbitalMultiplier::operator()(Orbital &Vphi, Potential &V, Orbital &phi) {
    if (Vphi.hasReal() or Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    Vphi.real = calcRealPart(V, phi);
    Vphi.imag = calcImagPart(V, phi, false);
}

void OrbitalMultiplier::adjoint(Orbital &Vphi, Potential &V, Orbital &phi) {
    if (Vphi.hasReal() or Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    Vphi.real = calcRealPart(V, phi);
    Vphi.imag = calcImagPart(V, phi, true);
}

// phi_ab = c * phi_a * phi_b
void OrbitalMultiplier::operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    phi_ab.real = calcRealPart(c, phi_a, phi_b);
    phi_ab.imag = calcImagPart(c, phi_a, phi_b, false);
}

// phi_ab = c * phi_a^dag * phi_b
void OrbitalMultiplier::adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    phi_ab.real = calcRealPart(c, phi_a, phi_b);
    phi_ab.imag = calcImagPart(c, phi_a, phi_b, true);
}

FunctionTree<3>* OrbitalMultiplier::calcRealPart(Potential &V, Orbital &phi) {
    FunctionTreeVector<3> vec;
    if (V.hasReal() and phi.hasReal()) {
        FunctionTree<3> *tree = this->grid(phi.re());
        this->mult(*tree, 1.0, V.re(), phi.re(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasImag()) {
        FunctionTree<3> *tree = this->grid(phi.im());
        this->mult(*tree, 1.0, V.im(), phi.im(), 0);
        vec.push_back(-1.0, tree);
    }
    FunctionTree<3> *real = 0;
    if (vec.size() == 1) {
        real = vec[0];
        vec.clear();
    }
    if (vec.size() == 2) {
        real = this->grid(vec);
        this->add(*real, vec, 0);
        vec.clear(true);
    }
    return real;
}

FunctionTree<3>* OrbitalMultiplier::calcImagPart(Potential &V, Orbital &phi, bool adjoint) {
    FunctionTreeVector<3> vec;
    if (V.hasReal() and phi.hasImag()) {
        FunctionTree<3> *tree = this->grid(phi.im());
        this->mult(*tree, 1.0, V.re(), phi.im(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasReal()) {
        FunctionTree<3> *tree = this->grid(phi.re());
        this->mult(*tree, 1.0, V.im(), phi.re(), 0);
        if (adjoint) {
            vec.push_back(-1.0, tree);
        } else {
            vec.push_back(1.0, tree);
        }
    }
    FunctionTree<3> *imag = 0;
    if (vec.size() == 1) {
        imag = vec[0];
        vec.clear();
    }
    if (vec.size() == 2) {
        imag = this->grid(vec);
        this->add(*imag, vec, 0);
        vec.clear(true);
    }
    return imag;
}
FunctionTree<3>* OrbitalMultiplier::calcRealPart(double c, Orbital &phi_a, Orbital &phi_b) {
    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasReal()) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, phi_a.re());
        this->grid(*tree, phi_b.re());
        this->mult(*tree, c, phi_a.re(), phi_b.re(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasImag()) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, phi_a.im());
        this->grid(*tree, phi_b.im());
        this->mult(*tree, c, phi_a.im(), phi_b.im(), 0);
        vec.push_back(-1.0, tree);
    }
    FunctionTree<3> *real = 0;
    if (vec.size() == 1) {
        real = vec[0];
        vec.clear();
    }
    if (vec.size() == 2) {
        real = this->grid(vec);
        this->add(*real, vec, 0);
        vec.clear(true);
    }
    return real;
}

FunctionTree<3>* OrbitalMultiplier::calcImagPart(double c, Orbital &phi_a, Orbital &phi_b, bool adjoint) {
    FunctionTreeVector<3> vec;
    if (phi_a.hasReal() and phi_b.hasImag()) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, phi_a.re());
        this->grid(*tree, phi_b.im());
        this->mult(*tree, c, phi_a.re(), phi_b.im(), 0);
        vec.push_back(1.0, tree);
    }
    if (phi_a.hasImag() and phi_b.hasReal()) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, phi_a.im());
        this->grid(*tree, phi_b.re());
        this->mult(*tree, c, phi_a.im(), phi_b.re(), 0);
        if (adjoint) {
            vec.push_back(-1.0, tree);
        } else {
            vec.push_back(1.0, tree);
        }
    }
    FunctionTree<3> *imag = 0;
    if (vec.size() == 1) {
        imag = vec[0];
        vec.clear();
    }
    if (vec.size() == 2) {
        imag = this->grid(vec);
        this->add(*imag, vec, 0);
        vec.clear(true);
    }
    return imag;
}
