#include "ExchangePotential.h"
#include "OrbitalAdder.h"
#include "OrbitalMultiplier.h"
#include "PoissonOperator.h"
#include "MWConvolution.h"

using namespace std;
using namespace Eigen;

ExchangePotential::ExchangePotential(PoissonOperator &P,
                                     OrbitalVector &phi,
                                     double x_fac)
        : ExchangeOperator(P, phi, x_fac),
          exchange(phi) {
}

void ExchangePotential::setup(double prec) {
    setApplyPrec(prec);
    calcInternalExchange();
}

void ExchangePotential::clear() {
    this->exchange.clear();
    clearApplyPrec();
}

void ExchangePotential::rotate(MatrixXd &U) {
    // negative prec means precomputed exchange is cleared and should therefore not be rotated
    if (this->apply_prec > 0.0) {
        OrbitalAdder add(this->apply_prec, this->max_scale);
        add.rotate(this->exchange, U);
    }
}

Orbital* ExchangePotential::operator() (Orbital &phi_p) {
    Orbital *preOrb = testPreComputed(phi_p);
    if (preOrb != 0) {
        println(1, "Precomputed exchange");
        return preOrb;
    }
    return calcExchange(phi_p);
}

Orbital* ExchangePotential::adjoint(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

/** Compute exchange on the fly */
Orbital* ExchangePotential::calcExchange(Orbital &phi_p) {
    Timer timer;

    OrbitalAdder add(this->apply_prec, this->max_scale);
    OrbitalMultiplier mult(this->apply_prec, this->max_scale);
    MWConvolution<3> apply(this->apply_prec, this->max_scale);
    PoissonOperator &P = *this->poisson;

    vector<complex<double> > coef_vec;
    vector<Orbital *> orb_vec;

    int nOrbs = this->orbitals->size();
    for (int i = 0; i < nOrbs; i++) {
        Orbital &phi_i = this->orbitals->getOrbital(i);

        Orbital *phi_ip = new Orbital(phi_p);
        mult.adjoint(*phi_ip, 1.0, phi_i, phi_p);

        Orbital *V_ip = new Orbital(phi_p);
        if (phi_ip->hasReal()) {
            V_ip->allocReal();
            apply(V_ip->real(), P, phi_ip->real());
        }
        if (phi_ip->hasImag()) {
            V_ip->allocImag();
            apply(V_ip->imag(), P, phi_ip->imag());
        }
        delete phi_ip;

        double spinFactor = phi_i.getExchangeFactor(phi_p);
        double fac = - spinFactor * (this->x_factor / phi_i.getSquareNorm());

        Orbital *phi_iip = new Orbital(phi_p);
        mult(*phi_iip, fac, phi_i, *V_ip);
        delete V_ip;

        coef_vec.push_back(1.0);
        orb_vec.push_back(phi_iip);
    }
    Orbital *ex_p = new Orbital(phi_p);
    add(*ex_p, coef_vec, orb_vec, true);

    for (int i = 0; i < orb_vec.size(); i++) delete orb_vec[i];
    orb_vec.clear();

    timer.stop();
    double n = ex_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied exchange", n, t);

    return ex_p;
}

/** Precompute the internal exchange */
void ExchangePotential::calcInternalExchange() {
    Timer timer;
    int nOrbs = this->orbitals->size();
    for (int i = 0; i < nOrbs; i++) {
        calcInternal(i);
        for (int j = 0; j < i; j++) {
            calcInternal(i,j);
        }
    }

    int n = 0;
    for (int i = 0; i < nOrbs; i++) {
        Orbital &ex_i = this->exchange.getOrbital(i);
        this->tot_norms(i) = sqrt(ex_i.getSquareNorm());
        n = max(n, ex_i.getNNodes());
    }

    timer.stop();
    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Hartree-Fock exchange", n, t);
}

void ExchangePotential::calcInternal(int i) {
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;
    Orbital &phi_i = this->orbitals->getOrbital(i);

    double prec = getScaledPrecision(i, i);
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ii = phi_i^dag * phi_i
    Orbital *phi_ii = new Orbital(phi_i);
    mult.adjoint(*phi_ii, 1.0, phi_i, phi_i);

    // compute V_ii = P[phi_ii]
    Orbital *V_ii = new Orbital(phi_i);
    if (phi_ii->hasReal()) {
        V_ii->allocReal();
        apply(V_ii->real(), P, phi_ii->real());
    }
    if (phi_ii->hasImag()) {
        V_ii->allocImag();
        apply(V_ii->imag(), P, phi_ii->imag());
    }
    if (phi_ii != 0) delete phi_ii;

    double fac_iii = -(this->x_factor/phi_i.getSquareNorm());

    // compute phi_iii = phi_i * V_ii
    Orbital *phi_iii = new Orbital(phi_i);
    mult(*phi_iii, fac_iii, phi_i, *V_ii);
    this->part_norms(i,i) = sqrt(phi_iii->getSquareNorm());
    if (V_ii != 0) delete V_ii;

    // compute x_i += phi_iii
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, 1.0, *phi_iii);
    if (phi_iii != 0) delete phi_iii;
}

void ExchangePotential::calcInternal(int i, int j) {
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;
    Orbital &phi_i = this->orbitals->getOrbital(i);
    Orbital &phi_j = this->orbitals->getOrbital(j);

    double i_factor = phi_i.getExchangeFactor(phi_j);
    double j_factor = phi_j.getExchangeFactor(phi_i);
    if (i_factor < MachineZero or j_factor < MachineZero) {
        this->part_norms(i,j) = 0.0;
        return;
    }

    double prec = getScaledPrecision(i, j);
    if (prec > 1.0e00) return;
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ij = phi_i^dag * phi_j
    Orbital *phi_ij = new Orbital(phi_i);
    mult.adjoint(*phi_ij, 1.0, phi_i, phi_j);

    // compute V_ij = P[phi_ij]
    Orbital *V_ij = new Orbital(phi_i);
    if (phi_ij->hasReal()) {
        V_ij->allocReal();
        apply(V_ij->real(), P, phi_ij->real());
    }
    if (phi_ij->hasImag()) {
        V_ij->allocImag();
        apply(V_ij->imag(), P, phi_ij->imag());
    }
    if (phi_ij != 0) delete phi_ij;

    // compute phi_jij = phi_j * V_ij
    double fac_jij = -(this->x_factor/phi_j.getSquareNorm());
    Orbital *phi_jij = new Orbital(phi_i);
    mult(*phi_jij, fac_jij, phi_j, *V_ij);
    this->part_norms(j,i) = sqrt(phi_jij->getSquareNorm());

    // compute phi_iij = phi_i * V_ij
    double fac_iij = -(this->x_factor/phi_i.getSquareNorm());
    Orbital *phi_iij = new Orbital(phi_i);
    mult(*phi_iij, fac_iij, phi_i, *V_ij);
    this->part_norms(i,j) = sqrt(phi_iij->getSquareNorm());

    if (V_ij != 0) delete V_ij;

    // compute x_i += phi_jij
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, i_factor, *phi_jij);
    if (phi_jij != 0) delete phi_jij;

    // compute x_j += phi_iij
    Orbital &ex_j = this->exchange.getOrbital(j);
    add.inPlace(ex_j, j_factor, *phi_iij);
    if (phi_iij != 0) delete phi_iij;
}

Orbital* ExchangePotential::testPreComputed(Orbital &phi_p) {
    int nOrbs = this->orbitals->size();
    for (int i = 0; i < nOrbs; i++) {
        Orbital &phi_i  = this->orbitals->getOrbital(i);
        Orbital *ex_i = this->exchange[i];
        if (&phi_i == &phi_p and ex_i != 0) {
            Orbital *result = new Orbital(phi_p);
            // Deep copy of orbital
            OrbitalAdder add(this->apply_prec, this->max_scale);
            add.inPlace(*result, 1.0, *ex_i);
            return result;
        }
    }
    return 0;
}
