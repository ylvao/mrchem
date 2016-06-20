#include "DensityProjector.h"
#include "OrbitalVector.h"
#include "Density.h"

using namespace std;

void DensityProjector::setPrecision(double prec) {
    this->clean.setPrecision(prec);
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void DensityProjector::operator()(Density &rho, Orbital &phi) {
    if (rho.total != 0) MSG_ERROR("Density not empty");
    if (rho.alpha != 0) MSG_ERROR("Density not empty");
    if (rho.beta != 0) MSG_ERROR("Density not empty");

    double occ = 1.0;
    if (not rho.spin) occ = (double) phi.getOccupancy();

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = this->mult(occ, phi.re(), phi.re());
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = this->mult(occ, phi.im(), phi.im());
        sum_vec.push_back(imag_2);
    }

    if (rho.spin) {
        if (phi.getSpin() == Paired) {
            rho.alpha = this->add(sum_vec);
            rho.beta = this->add(sum_vec);
        }
        if (phi.getSpin() == Alpha) {
            rho.alpha = this->add(sum_vec);
            rho.beta = this->grid(*rho.alpha);
            rho.beta->setZero();
        }
        if (phi.getSpin() == Beta) {
            rho.beta = this->add(sum_vec);
            rho.alpha = this->grid(*rho.beta);
            rho.alpha->setZero();
        }
        rho.total = this->add(1.0, *rho.alpha, 1.0, *rho.beta);
    } else {
        rho.total =  this->add(sum_vec);
        rho.alpha = 0;
        rho.beta = 0;
    }
    sum_vec.clear(true);
}

void DensityProjector::operator()(Density &rho, OrbitalVector &phi) {
    if (rho.total != 0) MSG_ERROR("Density not empty");
    if (rho.alpha != 0) MSG_ERROR("Density not empty");
    if (rho.beta != 0) MSG_ERROR("Density not empty");

    FunctionTreeVector<3> total_vec, alpha_vec, beta_vec;
    vector<Density *> dens_vec;
    for (int i = 0; i < phi.size(); i++) {
        Orbital &phi_i = phi.getOrbital(i);
        Density *rho_i = new Density(rho);
        (*this)(*rho_i, phi_i);
        dens_vec.push_back(rho_i);
        if (rho_i->total != 0) total_vec.push_back(rho_i->total);
        if (rho_i->alpha != 0) alpha_vec.push_back(rho_i->alpha);
        if (rho_i->beta != 0) beta_vec.push_back(rho_i->beta);
    }
    if (not rho.spin) {
        if (total_vec.size() > 0) rho.total = this->add(total_vec);
        rho.alpha = 0;
        rho.beta = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }

    for (int i = 0; i < dens_vec.size(); i++) {
        dens_vec[i]->clear();
        delete dens_vec[i];
        dens_vec[i] = 0;
    }
}
