#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

OrbitalAdder::OrbitalAdder(double pr)
    : add(*MRA, pr),
      grid(*MRA) {
}

void OrbitalAdder::operator()(Orbital &phi_ab,
                              double a, Orbital &phi_a,
                              double b, Orbital &phi_b,
                              bool union_grid) {
    double prec = this->add.getPrecision();
    if (not union_grid and prec < 0.0) MSG_ERROR("Nagative adaptive prec");
    if (phi_ab.hasReal() or phi_ab.hasImag()) MSG_ERROR("Orbital not empty");

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    if (phi_a.hasReal()) rvec.push_back(a, phi_a.real);
    if (phi_b.hasReal()) rvec.push_back(b, phi_b.real);

    if (phi_a.hasImag()) ivec.push_back(a, phi_a.imag);
    if (phi_b.hasImag()) ivec.push_back(b, phi_b.imag);

    if (rvec.size() > 0) {
        if (union_grid) {
            phi_ab.real = this->grid(rvec);
            this->add(*phi_ab.real, rvec, 0);
        } else {
            phi_ab.real = this->grid();
            this->add(*phi_ab.real, rvec);
        }
    }
    if (ivec.size() > 0) {
        if (union_grid) {
            phi_ab.imag = this->grid(ivec);
            this->add(*phi_ab.imag, ivec, 0);
        } else {
            phi_ab.imag = this->grid();
            this->add(*phi_ab.imag, ivec);
        }
    }
}

void OrbitalAdder::operator()(Orbital &out,
                              std::vector<double> &coefs,
                              std::vector<Orbital *> &orbs,
                              bool union_grid) {
    double prec = this->add.getPrecision();
    if (not union_grid and prec < 0.0) MSG_ERROR("Negative adaptive prec");
    if (out.hasReal() or out.hasImag()) MSG_ERROR("Orbital not empty");
    if (coefs.size() != orbs.size()) MSG_ERROR("Invalid arguments");

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;
    for (int i = 0; i < orbs.size(); i++) {
        if (orbs[i]->hasReal()) rvec.push_back(coefs[i], orbs[i]->real);
        if (orbs[i]->hasImag()) ivec.push_back(coefs[i], orbs[i]->imag);
    }

    if (rvec.size() > 0) {
        if (union_grid) {
            out.real = this->grid(rvec);
            this->add(*out.real, rvec, 0);
        } else {
            out.real = this->grid();
            this->add(*out.real, rvec);
        }
    }
    if (ivec.size() > 0) {
        if (union_grid) {
            out.imag = this->grid(ivec);
            this->add(*out.imag, ivec, 0);
        } else {
            out.imag = this->grid();
            this->add(*out.imag, ivec);
        }
    }
}

void OrbitalAdder::operator()(OrbitalVector &out,
                              double a, OrbitalVector &inp_a,
                              double b, OrbitalVector &inp_b,
                              bool union_grid) {
    if (out.size() != inp_a.size()) MSG_ERROR("Invalid arguments");
    if (out.size() != inp_b.size()) MSG_ERROR("Invalid arguments");

    for (int i = 0; i < out.size(); i++) {
        Orbital &out_i = out.getOrbital(i);
        Orbital &aInp_i = inp_a.getOrbital(i);
        Orbital &bInp_i = inp_b.getOrbital(i);
        (*this)(out_i, a, aInp_i, b, bInp_i, union_grid);
    }
}

void OrbitalAdder::operator()(Orbital &out,
                              const VectorXd &c,
                              OrbitalVector &inp,
                              bool union_grid) {
    double prec = this->add.getPrecision();
    if (not union_grid and prec < 0.0) MSG_ERROR("Nagative adaptive prec");
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

    if (rvec.size() > 0) {
        if (union_grid) {
            out.real = this->grid(rvec);
            this->add(*out.real, rvec, 0);
        } else {
            out.real = this->grid();
            this->add(*out.real, rvec);
        }
    }
    if (ivec.size() > 0) {
        if (union_grid) {
            out.imag = this->grid(ivec);
            this->add(*out.imag, ivec, 0);
        } else {
            out.imag = this->grid();
            this->add(*out.imag, ivec);
        }
    }
}

void OrbitalAdder::rotate(OrbitalVector &out, const MatrixXd &U, OrbitalVector &inp) {
    if (U.cols() != inp.size()) MSG_ERROR("Invalid arguments");
    if (U.rows() < out.size()) MSG_ERROR("Invalid arguments");

    for (int i = 0; i < out.size(); i++) {
        const VectorXd &c = U.row(i);
        Orbital &out_i = out.getOrbital(i);
        (*this)(out_i, c, inp, false); // Adaptive grids
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
    Orbital tmp(out);//shallow copy
    (*this)(tmp, 1.0, out, c, inp, true); // Union grid
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
