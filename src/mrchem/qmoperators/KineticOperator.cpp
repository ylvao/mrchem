#include "KineticOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

KineticOperator::KineticOperator(double build_prec,
                                 const MultiResolutionAnalysis<3> &mra)
        : QMOperator(mra),
          momentum_x(build_prec, mra, 0),
          momentum_y(build_prec, mra, 1),
          momentum_z(build_prec, mra, 2) {
}

KineticOperator::~KineticOperator() {
}

void KineticOperator::setup(double prec) {
    this->apply_prec = prec;
    this->momentum_x.setup(prec);
    this->momentum_y.setup(prec);
    this->momentum_z.setup(prec);
}

void KineticOperator::clear() {
    this->apply_prec = -1.0;
    this->momentum_x.clear();
    this->momentum_y.clear();
    this->momentum_z.clear();
}

Orbital* KineticOperator::operator() (Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

Orbital* KineticOperator::adjoint(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

double KineticOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    Orbital *xOrb_j = this->momentum_x(orb_j);
    Orbital *xOrb_i = this->momentum_x(orb_i);
    complex<double> T_x = xOrb_i->dot(*xOrb_j);
    delete xOrb_i;
    delete xOrb_j;

    Orbital *yOrb_j = this->momentum_y(orb_j);
    Orbital *yOrb_i = this->momentum_y(orb_i);
    complex<double> T_y = yOrb_i->dot(*yOrb_j);
    delete yOrb_i;
    delete yOrb_j;

    Orbital *zOrb_j = this->momentum_z(orb_j);
    Orbital *zOrb_i = this->momentum_z(orb_i);
    complex<double> T_z = zOrb_i->dot(*zOrb_j);
    delete zOrb_i;
    delete zOrb_j;

    complex<double> T_tot = T_x + T_y + T_z;
    if (T_tot.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return 0.5*T_tot.real();
}

double KineticOperator::adjoint(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd KineticOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    Timer timer;
    timer.restart();
    TelePrompter::printHeader(1, "Compute Kinetic Matrix Elements");

    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd T_x = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_y = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_z = MatrixXcd::Zero(Ni, Nj);
    {
        Timer timer;
        timer.restart();
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *xOrb_j = this->momentum_x(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *xOrb_i = this->momentum_x(orb_i);
                T_x(i,j) = xOrb_i->dot(*xOrb_j);
                delete xOrb_i;
            }
            delete xOrb_j;
        }
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_x", t);
    }
    {
        Timer timer;
        timer.restart();
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *yOrb_j = this->momentum_y(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *yOrb_i = this->momentum_y(orb_i);
                T_y(i,j) = yOrb_i->dot(*yOrb_j);
                delete yOrb_i;
            }
            delete yOrb_j;
        }
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_y", t);
    }
    {
        Timer timer;
        timer.restart();
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *zOrb_j = this->momentum_z(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *zOrb_i = this->momentum_z(orb_i);
                T_z(i,j) = zOrb_i->dot(*zOrb_j);
                delete zOrb_i;
            }
            delete zOrb_j;
        }
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_z", t);
    }
    MatrixXcd T_tot = T_x + T_y + T_z;

    TelePrompter::printFooter(1, timer, 2);
    if (T_tot.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return 0.5*T_tot.real();
}

MatrixXd KineticOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

int KineticOperator::printTreeSizes() const {
    println(0, " KineticOperator   " << setw(15) << 0 << setw(25) << 0);
    return 0;
}
