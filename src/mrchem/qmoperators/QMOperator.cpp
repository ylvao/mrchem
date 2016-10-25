#include "QMOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

QMOperator::QMOperator()
    : apply_prec(-1.0),
      clean(*MRA, -1.0),
      grid(*MRA){
}

void QMOperator::setup(double prec) {
    this->apply_prec = prec;
    this->clean.setPrecision(prec);
}
void QMOperator::clear() {
    this->apply_prec = -1.0;
    this->clean.setPrecision(-1.0);
}

/**

Calculates the matrix element of the operator between two orbitals:
\f$ \left\langle \phi_i|O|\phi_j \right\rangle\f$

The operator is applied "to the right" and then the inner product with the bra is taken.
 */
double QMOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    Orbital *operOrb = (*this)(orb_j);
    complex<double> result = orb_i.dot(*operOrb);
    delete operOrb;

    if (result.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return result.real();
}

/**

Calculates the matrix element of the adjoint operator between two orbitals:
\f$ \left\langle \phi_i|O^\dagger|\phi_j \right\rangle\f$
The operator is applied "to the left" and then the inner product with the ket is taken.

 */
double QMOperator::adjoint(Orbital &orb_i, Orbital &orb_j) {
    Orbital *operOrb = (*this).adjoint(orb_j);
    complex<double> result = orb_i.dot(*operOrb);
    delete operOrb;

    if (result.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return result.real();
}

/**

Calculates the matrix representation of the operator: \f$O\f$, given two sets of orbitals (bra and ket sides):
\f$ O^\dagger_{ij} = \left\langle \phi_i|O^\dagger|\phi_j \right\rangle\f$

 @param i_orbs bra functions
 @param j_orbs ket functions

 */
MatrixXd QMOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    TelePrompter::printHeader(1, "Compute Matrix Element");
    Timer timer;
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd result = MatrixXcd::Zero(Ni,Nj);
    for (int j = 0; j < Nj; j++) {
        Orbital &orb_j = j_orbs.getOrbital(j);
        Orbital *operOrb = (*this)(orb_j);
        for (int i = 0; i < Ni; i++) {
            Orbital &orb_i = i_orbs.getOrbital(i);
            result(i,j) = orb_i.dot(*operOrb);
        }
        delete operOrb;
    }
    timer.stop();
    TelePrompter::printFooter(1, timer, 2);
    if (result.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return result.real();
}

/**

Calculates the matrix representation of the adjoint operator: \f$O\f$, given two sets of orbitals (bra and ket sides):
\f$ O_{ij} = \left\langle \phi_i|O|\phi_j \right\rangle\f$

 @param i_orbs bra functions
 @param j_orbs ket functions

 */
MatrixXd QMOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    Timer timer;
    TelePrompter::printHeader(1, "Compute Matrix Element");
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd result = MatrixXcd::Zero(Ni,Nj);
    for (int i = 0; i < Ni; i++) {
        Orbital &orb_i = i_orbs.getOrbital(i);
        Orbital *operOrb = (*this)(orb_i);
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            result(i,j) = operOrb->dot(orb_j);
        }
        delete operOrb;
    }
    timer.stop();
    TelePrompter::printFooter(1, timer, 2);
    if (result.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return result.real();
}

/**
   Linear response property integrals: \f$ \eta_p \f$ is the occupation number, \f$\phi_p\f$ is the orbital \f$x_p\f$ and \f$y_p\f$ are the corresponding components of the response vectors.

   \f$ \eta_p ( \left\langle \phi_p|O^\dagger|x_p\right\rangle +
              \left\langle y_p |O|\phi_p \right\rangle ) \f$

 */
double QMOperator::calcProperty(Orbital &phi_p, Orbital *x_p, Orbital *y_p) {
    if (x_p == 0) x_p = &phi_p;
    if (y_p == 0) y_p = x_p;

    double occ = (double) phi_p.getOccupancy();
    double result_1 = (*this)(phi_p, *x_p);
    double result_2 = (*this)(*y_p, phi_p);
    return occ * (result_1 + result_2);
}

/**

   Linear response property: \f$ \eta_p \f$ is the occupation number, \f$\phi_p\f$ is the orbital \f$x_p\f$ and \f$y_p\f$ are the corresponding components of the response vectors.

   \f$ \sum_p \eta_p ( \left\langle \phi_p|O^\dagger |x_p\right\rangle +
              \left\langle y_p |O|\phi_p \right\rangle ) \f$

 */
double QMOperator::calcProperty(OrbitalVector &phi,
                                OrbitalVector *x,
                                OrbitalVector *y) {
    if (x == 0) x = &phi;
    if (y == 0) y = x;

    double result = 0.0;
    int nOrbs = phi.size();
    for (int i = 0; i < nOrbs; i++) {
        Orbital &phi_i = phi.getOrbital(i);
        Orbital &x_i = x->getOrbital(i);
        Orbital &y_i = y->getOrbital(i);
        result += calcProperty(phi_i, &x_i, &y_i);
    }
    return result;
}

int QMOperator::printTreeSizes() const {
    println(0, " QMOperator        " << setw(15) << 0 << setw(25) << 0);
    return 0;
}
