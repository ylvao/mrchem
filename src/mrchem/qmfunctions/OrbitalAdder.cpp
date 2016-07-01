#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"

using namespace std;
using namespace Eigen;

void OrbitalAdder::operator()(Orbital &phi_ab,
                              double a, Orbital &phi_a,
                              double b, Orbital &phi_b) {
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    if (phi_a.hasReal()) rvec.push_back(a, phi_a.real);
    if (phi_b.hasReal()) rvec.push_back(b, phi_b.real);

    if (phi_a.hasImag()) ivec.push_back(a, phi_a.imag);
    if (phi_b.hasImag()) ivec.push_back(b, phi_b.imag);

    // Fixed union grids
    if (rvec.size() > 0) {
        phi_ab.real = this->grid(rvec);
        this->add(*phi_ab.real, rvec, 0);
    }
    if (ivec.size() > 0) {
        phi_ab.imag = this->grid(ivec);
        this->add(*phi_ab.imag, ivec, 0);
    }
}

void OrbitalAdder::operator()(Orbital &out,
                              std::vector<double> &coefs,
                              std::vector<Orbital *> &orbs) {
    if (out.hasReal() or out.hasImag()) MSG_ERROR("Orbital not empty");
    if (coefs.size() != orbs.size()) MSG_ERROR("Invalid arguments");
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;
    for (int i = 0; i < orbs.size(); i++) {
        if (orbs[i]->hasReal()) rvec.push_back(coefs[i], orbs[i]->real);
        if (orbs[i]->hasImag()) ivec.push_back(coefs[i], orbs[i]->imag);
    }

    // Adaptive grids
    if (rvec.size() > 2) {
        out.real = this->grid();
        this->add(*out.real, rvec);
    } else if (rvec.size() > 0) {
        out.real = this->grid(rvec);
        this->add(*out.real, rvec, 0);
    }
    if (ivec.size() > 2) {
        out.imag = this->grid();
        this->add(*out.imag, ivec);
    } else if (ivec.size() > 0) {
        out.imag = this->grid(ivec);
        this->add(*out.imag, ivec, 0);
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

void OrbitalAdder::operator()(Orbital &out, const VectorXd &c, OrbitalVector &inp) {
    if (c.size() != inp.size()) MSG_ERROR("Invalid arguments");
    if (out.hasReal() or out.hasImag()) MSG_ERROR("Output not empty");

    double thrs = MachineZero;
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;
    for (int i = 0; i < inp.size(); i++) {
        double c_i = c(i);
        Orbital &phi_i = inp.getOrbital(i);
        if (phi_i.hasReal() and fabs(c_i) > thrs) rvec.push_back(c_i, phi_i.real);
        if (phi_i.hasImag() and fabs(c_i) > thrs) ivec.push_back(c_i, phi_i.imag);
    }

    // Adaptive grids
    if (rvec.size() > 2) {
        out.real = this->grid();
        this->add(*out.real, rvec);
    } else if (rvec.size() > 0) {
        out.real = this->grid(rvec);
        this->add(*out.real, rvec, 0);
    }
    if (ivec.size() > 2) {
        out.imag = this->grid();
        this->add(*out.imag, ivec);
    } else if (ivec.size() > 0) {
        out.imag = this->grid(ivec);
        this->add(*out.imag, ivec, 0);
    }
}

void OrbitalAdder::rotate(OrbitalVector &out, const MatrixXd &U, OrbitalVector &inp) {
    if (out.size() != inp.size()) MSG_ERROR("Invalid arguments");
    if (out.size() != U.rows()) MSG_ERROR("Invalid arguments");
    for (int i = 0; i < out.size(); i++) {
        const VectorXd &c = U.row(i);
        Orbital &out_i = out.getOrbital(i);
        (*this)(out_i, c, inp);
    }
}

/** In place rotation of orbital vector */
void OrbitalAdder::rotate(OrbitalVector &out, const MatrixXd &U) {
    OrbitalVector tmp(out);
    rotate(tmp, U, out);
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

void OrbitalAdder::inPlace(Orbital &out, double c, Orbital &inp) {
    Orbital tmp(out);
    (*this)(tmp, 1.0, out, c, inp);
    out.clear();
    out.real = tmp.real;
    out.imag = tmp.imag;
    tmp.real = 0;
    tmp.imag = 0;
}

void OrbitalAdder::inPlace(OrbitalVector &out, double c, OrbitalVector &inp) {
    if (out.size() != inp.size()) MSG_ERROR("Invalid arguments");

    for (int i = 0; i < out.size(); i++) {
        Orbital &out_i = out.getOrbital(i);
        Orbital &inp_i = inp.getOrbital(i);
        this->inPlace(out_i, c, inp_i);
    }
}
