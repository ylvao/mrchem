#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"

using namespace std;
using namespace Eigen;

void OrbitalAdder::operator()(Orbital &phi_ab,
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

void OrbitalAdder::operator()(Orbital &out,
                              std::vector<double> &coefs,
                              std::vector<Orbital *> &orbs) {
    if (out.hasReal() or out.hasImag()) MSG_ERROR("Orbital not empty");
    if (coefs.size() != orbs.size()) MSG_ERROR("Invalid arguments");
    {
        FunctionTreeVector<3> vec;
        for (int i = 0; i < orbs.size(); i++) {
            if (orbs[i]->hasReal()) vec.push_back(coefs[i], orbs[i]->real);
        }
        out.real = (*this)(vec);
    }
    {
        FunctionTreeVector<3> vec;
        for (int i = 0; i < orbs.size(); i++) {
            if (orbs[i]->hasImag()) vec.push_back(coefs[i], orbs[i]->imag);
        }
        out.imag = (*this)(vec);
    }
}

void OrbitalAdder::operator()(OrbitalVector &out,
                              double a, OrbitalVector &inp_a,
                              double b, OrbitalVector &inp_b) {
    if (out.size() != inp_a.size()) MSG_ERROR("Invalid arguments");
    if (out.size() != inp_b.size()) MSG_ERROR("Invalid arguments");

    for (int i = 0; i < out.size(); i++) {
        Orbital &out_i = out.getOrbital(i);
        Orbital &aInp_i = inp_a.getOrbital(i);
        Orbital &bInp_i = inp_b.getOrbital(i);
        (*this)(out_i, a, aInp_i, b, bInp_i);
    }
}

void OrbitalAdder::inPlace(OrbitalVector &out, double c, OrbitalVector &inp) {
    if (out.size() != inp.size()) MSG_ERROR("Invalid arguments");

    for (int i = 0; i < out.size(); i++) {
        Orbital &out_i = out.getOrbital(i);
        Orbital &inp_i = inp.getOrbital(i);
        this->inPlace(out_i, c, inp_i);
    }
}

void OrbitalAdder::inPlace(Orbital &out, double c, Orbital &inp) {
    Orbital *tmp = new Orbital(out);
    (*this)(*tmp, 1.0, out, c, inp);
    out.clear();
    out.real = tmp->real;
    out.imag = tmp->imag;
    tmp->real = 0;
    tmp->imag = 0;
}


void OrbitalAdder::rotate(OrbitalVector &out, MatrixXd &U, OrbitalVector &inp) {
    NOT_IMPLEMENTED_ABORT;
}

void OrbitalAdder::operator()(Orbital &out, VectorXd &c, OrbitalVector &inp) {
    NOT_IMPLEMENTED_ABORT;
}
