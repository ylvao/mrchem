#include "OrbitalMultiplier.h"
#include "Orbital.h"
#include "Potential.h"

OrbitalMultiplier::OrbitalMultiplier(const MultiResolutionAnalysis<3> &mra, double pr)
    : add(mra, pr),
      mult(mra, pr),
      grid(mra) {
}

void OrbitalMultiplier::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void OrbitalMultiplier::operator()(Orbital &Vphi, double c, Potential &V, Orbital &phi) {
    if (Vphi.hasReal() or Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    Vphi.real = calcRealPart(c, V.real, V.imag, phi.real, phi.imag);
    Vphi.imag = calcImagPart(c, V.real, V.imag, phi.real, phi.imag, false);
}

void OrbitalMultiplier::adjoint(Orbital &Vphi, double c, Potential &V, Orbital &phi) {
    if (Vphi.hasReal() or Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    Vphi.real = calcRealPart(c, V.real, V.imag, phi.real, phi.imag);
    Vphi.imag = calcImagPart(c, V.real, V.imag, phi.real, phi.imag, true);
}

// phi_ab = c * phi_a * phi_b
void OrbitalMultiplier::operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    phi_ab.real = calcRealPart(c, phi_a.real, phi_a.imag, phi_b.real, phi_b.imag);
    phi_ab.imag = calcImagPart(c, phi_a.real, phi_a.imag, phi_b.real, phi_b.imag, false);
}

// phi_ab = c * phi_a^dag * phi_b
void OrbitalMultiplier::adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    phi_ab.real = calcRealPart(c, phi_a.real, phi_a.imag, phi_b.real, phi_b.imag);
    phi_ab.imag = calcImagPart(c, phi_a.real, phi_a.imag, phi_b.real, phi_b.imag, true);
}

FunctionTree<3>* OrbitalMultiplier::calcRealPart(double c,
                                                 FunctionTree<3> *re_a,
                                                 FunctionTree<3> *im_a,
                                                 FunctionTree<3> *re_b,
                                                 FunctionTree<3> *im_b) {
    FunctionTreeVector<3> vec;
    if (re_a != 0 and re_b != 0) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, *re_a);
        this->grid(*tree, *re_b);
        this->mult(*tree, c, *re_a, *re_b, 0);
        vec.push_back(1.0, tree);
    }
    if (im_a != 0 and im_b != 0) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, *im_a);
        this->grid(*tree, *im_b);
        this->mult(*tree, c, *im_a, *im_b, 0);
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

FunctionTree<3>* OrbitalMultiplier::calcImagPart(double c,
                                                 FunctionTree<3> *re_a,
                                                 FunctionTree<3> *im_a,
                                                 FunctionTree<3> *re_b,
                                                 FunctionTree<3> *im_b,
                                                 bool adjoint) {
    FunctionTreeVector<3> vec;
    if (re_a != 0 and im_b != 0) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, *re_a);
        this->grid(*tree, *im_b);
        this->mult(*tree, c, *re_a, *im_b, 0);
        vec.push_back(1.0, tree);
    }
    if (im_a != 0 and re_b != 0) {
        FunctionTree<3> *tree = this->grid();
        this->grid(*tree, *im_a);
        this->grid(*tree, *re_b);
        this->mult(*tree, c, *im_a, *re_b, 0);
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
