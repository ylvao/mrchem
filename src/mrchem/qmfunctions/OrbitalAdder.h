#ifndef ORBITALADDER_H
#define ORBITALADDER_H

#include "MWAdder.h"
#include "Orbital.h"

class OrbitalAdder : public MWAdder<3> {
public:
    OrbitalAdder(const MultiResolutionAnalysis<3> &mra, double pr = -1.0)
        : MWAdder(mra, pr) { }
    virtual ~OrbitalAdder() { }

    void operator()(Orbital &phi_ab,
                    double a, Orbital &phi_a,
                    double b, Orbital &phi_b) {
        if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
        {
            FunctionTreeVector<3> vec;
            if (phi_a.hasReal()) vec.push_back(a, phi_a.real);
            if (phi_b.hasReal()) vec.push_back(b, phi_b.real);
            phi_ab.real = (*this)(vec);
        }
        {
            FunctionTreeVector<3> vec;
            if (phi_a.hasImag()) vec.push_back(a, phi_a.imag);
            if (phi_b.hasImag()) vec.push_back(b, phi_b.imag);
            phi_ab.imag = (*this)(vec);
        }
    }

    using MWAdder<3>::operator();
};

#endif // ORBITALADDER_H
