#include "AngularOperator.h"
#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

AngularOperator::AngularOperator(int dir, double c, const double *o)
    : QMOperator(MRA->getMaxScale()),
      apply_dir(dir),
      coef(c),
      r_x(0, o),
      r_y(1, o),
      r_z(2, o),
      p_x(0),
      p_y(1),
      p_z(2) {
}

void AngularOperator::setup(double prec) {
    QMOperator::setup(prec);
    this->r_x.setup(prec);
    this->r_y.setup(prec);
    this->r_z.setup(prec);
    this->p_x.setup(prec);
    this->p_y.setup(prec);
    this->p_z.setup(prec);
}

void AngularOperator::clear() {
    this->r_x.clear();
    this->r_y.clear();
    this->r_z.clear();
    this->p_x.clear();
    this->p_y.clear();
    this->p_z.clear();
    QMOperator::clear();
}

Orbital* AngularOperator::operator()(Orbital &orb_p) {
    Timer timer;
    OrbitalAdder add(this->apply_prec);
    Orbital *result = new Orbital(orb_p);
    if (this->apply_dir == 0) {
        Timer t1;
        Orbital *dy = this->p_y(orb_p);
        Orbital *rzdy = this->r_z(*dy);
        delete dy;
        t1.stop();
        TelePrompter::printTree(2, "Mult rzdy", rzdy->getNNodes(), t1.getWallTime());

        Timer t2;
        Orbital *dz = this->p_z(orb_p);
        Orbital *rydz = this->r_y(*dz);
        delete dz;
        t2.stop();
        TelePrompter::printTree(2, "Mult rydz", rydz->getNNodes(), t2.getWallTime());

        add(*result, -this->coef, *rydz, this->coef, *rzdy, true);

        delete rzdy;
        delete rydz;
    } else if (this->apply_dir == 1) {
        Timer t1;
        Orbital *dx = this->p_x(orb_p);
        Orbital *rzdx = this->r_z(*dx);
        delete dx;
        t1.stop();
        TelePrompter::printTree(2, "Mult rzdx", rzdx->getNNodes(), t1.getWallTime());

        Timer t2;
        Orbital *dz = this->p_z(orb_p);
        Orbital *rxdz = this->r_x(*dz);
        delete dz;
        t2.stop();
        TelePrompter::printTree(2, "Mult rxdz", rxdz->getNNodes(), t2.getWallTime());

        add(*result, this->coef, *rxdz, -this->coef, *rzdx, true);

        delete rzdx;
        delete rxdz;
    } else if (this->apply_dir == 2) {
        Timer t1;
        Orbital *dx = this->p_x(orb_p);
        Orbital *rydx = this->r_y(*dx);
        delete dx;
        t1.stop();
        TelePrompter::printTree(2, "Mult rydx", rydx->getNNodes(), t1.getWallTime());

        Timer t2;
        Orbital *dy = this->p_y(orb_p);
        Orbital *rxdy = this->r_x(*dy);
        delete dy;
        t2.stop();
        TelePrompter::printTree(2, "Mult rxdy", rxdy->getNNodes(), t2.getWallTime());

        add(*result, -this->coef, *rxdy, this->coef, *rydx, true);

        delete rydx;
        delete rxdy;
    } else {
        MSG_ERROR("Invalid dir");
    }
    timer.stop();
    TelePrompter::printTree(1, "Angular operator", result->getNNodes(), timer.getWallTime());
    return result;
}

Orbital* AngularOperator::adjoint(Orbital &orb_p) {
    Orbital *result = (*this)(orb_p);
    return result;
}
