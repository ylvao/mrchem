#ifndef ORBITALMULTIPLIER_H
#define ORBITALMULTIPLIER_H

#include "MWMultiplier.h"
#include "MWAdder.h"
#include "Orbital.h"

class OrbitalMultiplier : public MWMultiplier<3> {
public:
    OrbitalMultiplier(const MultiResolutionAnalysis<3> &mra, double pr = -1.0)
        : MWMultiplier(mra, pr),
          add(mra) { }
    virtual ~OrbitalMultiplier() { }

    // phi_ab = c * phi_a * phi_b
    void operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
        if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
        phi_ab.real = calcRealPart(c, phi_a, phi_b);
        phi_ab.imag = calcImagPart(c, phi_a, phi_b, false);
    }

    // phi_ab = c * phi_a^dag * phi_b
    void adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b) {
        if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
        phi_ab.real = calcRealPart(c, phi_a, phi_b);
        phi_ab.imag = calcImagPart(c, phi_a, phi_b, true);
    }

protected:
    MWAdder<3> add;

    FunctionTree<3> *calcRealPart(double c, Orbital &phi_a, Orbital &phi_b) {
        FunctionTreeVector<3> vec;
        if (phi_a.hasReal() and phi_b.hasReal()) {
            FunctionTree<3> *tree = (*this)(c, phi_a.re(), phi_b.re());
            vec.push_back(1.0, tree);
        }
        if (phi_a.hasImag() and phi_b.hasImag()) {
            FunctionTree<3> *tree = (*this)(c, phi_a.im(), phi_b.im());
            vec.push_back(-1.0, tree);
        }
        FunctionTree<3> *real = 0;
        if (vec.size() == 1) {
            real = vec[0];
            vec.clear();
        }
        if (vec.size() == 2) {
            real = this->add(vec);
            vec.clear(true);
        }
        return real;
    }

    FunctionTree<3> *calcImagPart(double c, Orbital &phi_a, Orbital &phi_b, bool adjoint) {
        FunctionTreeVector<3> vec;
        if (phi_a.hasReal() and phi_b.hasImag()) {
            FunctionTree<3> *tree = (*this)(c, phi_a.re(), phi_b.im());
            vec.push_back(1.0, tree);
        }
        if (phi_a.hasImag() and phi_b.hasReal()) {
            FunctionTree<3> *tree = (*this)(c, phi_a.im(), phi_b.re());
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
            imag = this->add(vec);
            vec.clear(true);
        }
        return imag;
    }

    using MWMultiplier<3>::operator();
};


#endif // ORBITALMULTIPLIER_H
