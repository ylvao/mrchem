#include <fstream>

#include "Potential.h"
#include "FunctionTreeVector.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "MWAdder.h"
#include "MWMultiplier.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

Potential::Potential(const MultiResolutionAnalysis<3> &mra)
    : add(mra),
      mult(mra),
      real(0),
      imag(0) {
}

Potential::~Potential() {
    if (this->real != 0) MSG_ERROR("Operator not properly deallocated");
    if (this->imag != 0) MSG_ERROR("Operator not properly deallocated");
}

void Potential::setup(double prec) {
    this->apply_prec = prec;
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void Potential::clear() {
    this->apply_prec = -1.0;
    if (this->real != 0) delete this->real;
    if (this->imag != 0) delete this->imag;
    this->real = 0;
    this->imag = 0;
}


Orbital* Potential::operator() (Orbital &orb_p) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Timer timer;
    timer.restart();

    FunctionTree<3> *real_1 = 0;
    FunctionTree<3> *real_2 = 0;
    FunctionTree<3> *imag_1 = 0;
    FunctionTree<3> *imag_2 = 0;

    FunctionTreeVector<3> real_vec;
    FunctionTreeVector<3> imag_vec;

    Orbital *Vorb = new Orbital(orb_p);
    if (orb_p.real != 0) {
        if (this->real != 0) {
            FunctionTreeVector<3> tree_vec;
            tree_vec.push_back(this->real);
            tree_vec.push_back(orb_p.real);
            real_1 = this->mult(tree_vec);
            real_vec.push_back(1.0, real_1);
        }
        if (this->imag != 0) {
            FunctionTreeVector<3> tree_vec;
            tree_vec.push_back(this->imag);
            tree_vec.push_back(orb_p.real);
            imag_1 = this->mult(tree_vec);
            imag_vec.push_back(1.0, imag_1);
        }
    }
    if (orb_p.imag != 0) {
        if (this->real != 0) {
            FunctionTreeVector<3> tree_vec;
            tree_vec.push_back(this->real);
            tree_vec.push_back(orb_p.imag);
            imag_2 = this->mult(tree_vec);
            imag_vec.push_back(1.0, imag_2);
        }
        if (this->imag != 0) {
            FunctionTreeVector<3> tree_vec;
            tree_vec.push_back(this->imag);
            tree_vec.push_back(orb_p.imag);
            real_2 = this->mult(tree_vec);
            real_vec.push_back(-1.0, real_2);
        }
    }
    if (real_vec.size() > 0) {
        Vorb->real = this->add(real_vec);
    }
    if (imag_vec.size() > 0) {
        Vorb->imag = this->add(imag_vec);
    }

    if (real_1 != 0) delete real_1;
    if (real_2 != 0) delete real_2;
    if (imag_1 != 0) delete imag_1;
    if (imag_2 != 0) delete imag_2;

    timer.stop();
    int n = Vorb->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied potential", n, t);
    return Vorb;
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
