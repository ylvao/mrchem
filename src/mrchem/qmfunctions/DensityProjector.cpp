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
        FunctionTree<3> *real_2 = this->grid(phi.re());
        this->mult(*real_2, occ, phi.re(), phi.re(), 1);
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = this->grid(phi.im());
        this->mult(*imag_2, occ, phi.im(), phi.im(), 1);
        sum_vec.push_back(imag_2);
    }

    if (rho.spin) {
        if (phi.getSpin() == Paired) {
            rho.alpha = this->grid(sum_vec);
            rho.beta = this->grid(sum_vec);
            this->add(*rho.alpha, sum_vec, 0);
            this->add(*rho.beta, sum_vec, 0);
        }
        if (phi.getSpin() == Alpha) {
            rho.alpha = this->grid(sum_vec);
            this->add(*rho.alpha, sum_vec);
            rho.beta = this->grid(sum_vec);
            rho.beta->setZero();
        }
        if (phi.getSpin() == Beta) {
            rho.beta = this->grid(sum_vec);
            this->add(*rho.beta, sum_vec, 0);
            rho.alpha = this->grid(sum_vec);
            rho.alpha->setZero();
        }
        FunctionTreeVector<3> tot_vec;
        tot_vec.push_back(1.0, rho.alpha);
        tot_vec.push_back(1.0, rho.beta);
        rho.total = this->grid(tot_vec);
        this->add(*rho.total, tot_vec, 0);
    } else {
        rho.total = this->grid(sum_vec);
        this->add(*rho.total, sum_vec, 0);
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
        if (total_vec.size() > 0) {
            rho.total = this->grid(total_vec);
            this->add(*rho.total, total_vec, 0);
        }
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
