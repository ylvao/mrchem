#include "KineticOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

double KineticOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    Orbital *xOrb_j = p_x(orb_j);
    Orbital *xOrb_i = p_x(orb_i);
    complex<double> T_x = xOrb_i->dot(*xOrb_j);
    delete xOrb_i;
    delete xOrb_j;

    Orbital *yOrb_j = p_y(orb_j);
    Orbital *yOrb_i = p_y(orb_i);
    complex<double> T_y = yOrb_i->dot(*yOrb_j);
    delete yOrb_i;
    delete yOrb_j;

    Orbital *zOrb_j = p_z(orb_j);
    Orbital *zOrb_i = p_z(orb_i);
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
    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    Timer timer;
    TelePrompter::printHeader(1, "Compute Kinetic Matrix Elements");

    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd T_x = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_y = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_z = MatrixXcd::Zero(Ni, Nj);
    {
        Timer timer;
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *xOrb_j = p_x(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *xOrb_i = p_x(orb_i);
                T_x(i,j) = xOrb_i->dot(*xOrb_j);
                delete xOrb_i;
            }
            delete xOrb_j;
        }
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_x", t);
    }
    {
        Timer timer;
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *yOrb_j = p_y(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *yOrb_i = p_y(orb_i);
                T_y(i,j) = yOrb_i->dot(*yOrb_j);
                delete yOrb_i;
            }
            delete yOrb_j;
        }
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_y", t);
    }
    {
        Timer timer;
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *zOrb_j = p_z(orb_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &orb_i = i_orbs.getOrbital(i);
                Orbital *zOrb_i = p_z(orb_i);
                T_z(i,j) = zOrb_i->dot(*zOrb_j);
                delete zOrb_i;
            }
            delete zOrb_j;
        }
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_z", t);
    }
    MatrixXcd T_tot = T_x + T_y + T_z;

    timer.stop();
    TelePrompter::printFooter(1, timer, 2);
    if (T_tot.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return 0.5*T_tot.real();
}

MatrixXd KineticOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

