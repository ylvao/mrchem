#include "CoulombPotential.h"
#include "DensityProjector.h"
#include "PoissonOperator.h"
#include "MWConvolution.h"

using namespace std;

void CoulombPotential::setup(double prec) {
    QMOperator::setup(prec);
    this->potential.setup(prec);

    if (this->orbitals == 0) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &phi = *this->orbitals;
    Density &rho = this->density;
    QMPotential &V = this->potential;
    PoissonOperator &P = *this->poisson;

    MWConvolution<3> apply(this->apply_prec, this->max_scale);
    DensityProjector project(this->apply_prec, this->max_scale);

    Timer timer1;
    project(rho, phi);
    timer1.stop();
    double t1 = timer1.getWallTime();
    int n1 = rho.getNNodes();
    TelePrompter::printTree(0, "Coulomb density", n1, t1);

    Timer timer2;
    V.allocReal();
    apply(V.real(), P, rho.total());
    timer2.stop();
    int n2 = V.getNNodes();
    double t2 = timer2.getWallTime();
    TelePrompter::printTree(0, "Coulomb potential", n2, t2);
}

void CoulombPotential::clear() {
    this->density.clear();
    this->potential.clear();
    QMOperator::clear();
}

Orbital* CoulombPotential::operator() (Orbital &phi) {
    QMPotential &V = this->potential;
    return V(phi);
}

Orbital* CoulombPotential::adjoint(Orbital &phi) {
    QMPotential &V = this->potential;
    return V.adjoint(phi);
}
