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

void OrbitalAdder::operator()(Orbital &out, VectorXd &c, OrbitalVector &inp) {
    if (c.size() != inp.size()) MSG_ERROR("Invalid arguments");
    if (out.hasReal() or out.hasImag()) MSG_ERROR("Output not empty");

    FunctionTreeVector<3> real_vec;
    FunctionTreeVector<3> imag_vec;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &phi_i = inp.getOrbital(i);
        if (phi_i.hasReal()) real_vec.push_back(c(i), phi_i.real);
        if (phi_i.hasImag()) imag_vec.push_back(c(i), phi_i.imag);
    }
    if (real_vec.size() != 0) out.real = (*this)(real_vec);
    if (imag_vec.size() != 0) out.imag = (*this)(imag_vec);
}

void OrbitalAdder::rotate(OrbitalVector &out, MatrixXd &U, OrbitalVector &inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** In place rotation of orbital vector */
void OrbitalAdder::rotate(OrbitalVector &out, MatrixXd &U) {
    OrbitalVector tmp(out);
    for (int i = 0; i < out.size(); i++) {
        VectorXd c = U.row(i);
        Orbital &tmp_phi = tmp.getOrbital(i);
        (*this)(tmp_phi, c, out);
    }
    out.clear();
    for (int i = 0; i < out.size(); i++) {
        Orbital &tmp_phi = tmp.getOrbital(i);
        Orbital &out_phi = out.getOrbital(i);
        if (tmp_phi.hasReal()) {
            out_phi.real = tmp_phi.real;
            tmp_phi.real = 0;
        }
        if (tmp_phi.hasImag()) {
            out_phi.imag = tmp_phi.imag;
            tmp_phi.imag = 0;
        }
    }
    tmp.clear();
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
