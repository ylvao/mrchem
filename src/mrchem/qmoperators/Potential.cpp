#include <fstream>

#include "Potential.h"
#include "FunctionTreeVector.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

Potential::Potential(const MultiResolutionAnalysis<3> &mra)
    : QMOperator(mra),
      mult(mra, -1.0),
      real(0),
      imag(0) {
}

Potential::~Potential() {
    if (this->real != 0) MSG_ERROR("Operator not properly deallocated");
    if (this->imag != 0) MSG_ERROR("Operator not properly deallocated");
}

void Potential::setup(double prec) {
    QMOperator::setup(prec);
    this->mult.setPrecision(prec);
}

void Potential::clear() {
    if (this->real != 0) delete this->real;
    if (this->imag != 0) delete this->imag;
    this->real = 0;
    this->imag = 0;
    this->mult.setPrecision(-1.0);
    QMOperator::clear();
}


Orbital* Potential::operator() (Orbital &phi_p) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Timer timer;
    timer.restart();

    Potential &V = *this;
    Orbital *Vphi_p = new Orbital(phi_p);
    this->mult(*Vphi_p, 1.0, V, phi_p);

    timer.stop();
    int n = Vphi_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied potential", n, t);
    return Vphi_p;
}

Orbital* Potential::adjoint(Orbital &orb) {
    NOT_IMPLEMENTED_ABORT;
}

double Potential::operator() (Orbital &orb_i, Orbital &orb_j) {
    Orbital *operOrb = (*this)(orb_j);
    complex<double> result = orb_i.dot(*operOrb);
    if (result.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    delete operOrb;
    return result.real();
}

double Potential::adjoint(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd Potential::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    Timer timer;
    timer.restart();
    TelePrompter::printHeader(1, "Compute Potential Matrix Elements");
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd M = MatrixXcd::Zero(Ni,Nj);
    for (int j = 0; j < Nj; j++) {
        Orbital &orb_j = j_orbs.getOrbital(j);
        Orbital *operOrb = (*this)(orb_j);
        for (int i = 0; i < Ni; i++) {
            Orbital &orb_i = i_orbs.getOrbital(i);
            M(i,j) = orb_i.dot(*operOrb);
        }
        delete operOrb;
    }
    TelePrompter::printFooter(1, timer, 2);
    if (M.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return M.real();
}

MatrixXd Potential::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

int Potential::getNNodes() const {
    int nNodes = 0;
    if (this->real != 0) nNodes += this->real->getNNodes();
    if (this->imag != 0) nNodes += this->imag->getNNodes();
    return nNodes;
}

int Potential::printTreeSizes() const {
    int nNodes = this->getNNodes();
    println(0, " Potential         " << setw(15) << 1 << setw(25) << nNodes);
    return nNodes;
}
