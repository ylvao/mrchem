#include "ExchangePotential.h"

using namespace std;
using namespace Eigen;

ExchangePotential::ExchangePotential(double build_prec,
                                     OrbitalVector &phi,
                                     double x_fac)
        : ExchangeOperator(build_prec, phi, x_fac),
          exchange_0(phi) {
}

ExchangePotential::~ExchangePotential() {
}

void ExchangePotential::setup(double prec) {
    ExchangeOperator::setup(prec);
//    calcInternalExchange();
}

void ExchangePotential::clear() {
    this->exchange_0.clear();
    ExchangeOperator::clear();
}

void ExchangePotential::rotate(MatrixXd &U) {
    this->add.rotate(exchange_0, U);
}

Orbital* ExchangePotential::operator() (Orbital &phi_p) {
//    Orbital *preOrb = testPreComputed(phi_p);
//    if (preOrb != 0) {
//        println(1, "Precomputed exchange");
//        return preOrb;
//    }
    return calcExchange(phi_p);
}

Orbital* ExchangePotential::adjoint(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
//    return (*this)(phi_p);
}

/** Compute exchange on the fly */
Orbital* ExchangePotential::calcExchange(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
//    Timer timer;

//    int maxNodes = 0;
//    FunctionTreeVector<3> real_vec;
//    FunctionTreeVector<3> imag_vec;
//    int nOrbs = this->orbitals_0->size();
//    for (int i = 0; i < nOrbs; i++) {
//        Orbital &phi_i = this->orbitals_0->getOrbital(i);

//        Orbital *phi_ip = new Orbital(phi_p);
//        { // compute phi_ip = phi_i^dag * phi_p
//            this->mult.setPrecision(this->apply_prec);
//            this->mult.adjoint(*phi_ip, 1.0, phi_i, phi_p);
//            int nNodes = phi_ip->getNNodes();
//            double time = timer.getWallTime();
//            maxNodes = max(maxNodes, nNodes);
//            TelePrompter::printTree(2, "Multiply", nNodes, time);
//        }

//        Orbital *V_ip = new Orbital(phi_p);
//        { // compute V_ip = P[phi_ip]
//            Timer timer;
//            this->poisson.setPrecision(this->apply_prec);
//            if (phi_ip->hasReal()) {
//                V_ip->real = this->grid();
//                this->poisson(*V_ip->real, *phi_ip->real);
//            }
//            if (phi_ip->hasImag()) {
//                V_ip->imag = this->grid();
//                this->poisson(*V_ip->imag, *phi_ip->imag);
//            }
//            int nNodes = V_ip->getNNodes();
//            double time = timer.getWallTime();
//            maxNodes = max(maxNodes, nNodes);
//            TelePrompter::printTree(2, "Applying Poisson", nNodes, time);
//        }
//        delete phi_ip;

//        Orbital *phi_iip = new Orbital(phi_p);
//        { // compute phi_iij = phi_i * V_ip
//            Timer timer;
//            double fac = -(this->x_factor/phi_i.getSquareNorm());
//            this->mult(*phi_iip, fac, phi_i, *V_ip);
//            int nNodes = phi_iip->getNNodes();
//            double time = timer.getWallTime();
//            maxNodes = max(maxNodes, nNodes);
//            TelePrompter::printTree(2, "Multiply", nNodes, time);
//        }
//        delete V_ip;

//        if (phi_iip->hasReal()) real_vec.push_back(phi_iip->real);
//        if (phi_iip->hasImag()) imag_vec.push_back(phi_iip->imag);
//        phi_iip->real = 0;
//        phi_iip->imag = 0;
//        delete phi_iip;
//    }

//    Orbital *ex_p = new Orbital(phi_p);
//    if (real_vec.size() > 0) {
//        Timer timer;
//        ex_p->real = this->add(real_vec);
//        double time = timer.getWallTime();
//        double nNodes = ex_p->real->getNNodes();
//        TelePrompter::printTree(2, "Sum exchange", nNodes, time);
//    }
//    real_vec.clear(true);

//    if (imag_vec.size() > 0) {
//        Timer timer;
//        ex_p->imag = this->add(imag_vec);
//        double time = timer.getWallTime();
//        double nNodes = ex_p->imag->getNNodes();
//        TelePrompter::printTree(2, "Sum exchange", nNodes, time);
//    }
//    imag_vec.clear(true);

//    double nNodes = ex_p->getNNodes();
//    double time = timer.getWallTime();
//    TelePrompter::printTree(1, "Applied exchange", nNodes, time);
//    return ex_p;
}

/** Precompute the internal exchange */
void ExchangePotential::calcInternalExchange() {
    NOT_IMPLEMENTED_ABORT;
    //Timer timer;

    //int nNodes = 0;
    //int maxNodes = 0;

    //int nOrbs = this->orbitals_0->size();
    //for (int i = 0; i < nOrbs; i++) {
        //nNodes = calcInternal(i);
        //maxNodes = max(maxNodes, nNodes);
        //for (int j = 0; j < i; j++) {
            //nNodes = calcInternal(i,j);
            //maxNodes = max(maxNodes, nNodes);
        //}
    //}

    //for (int i = 0; i < nOrbs; i++) {
        //Orbital &ex_i = this->exchange_0.getOrbital(i);
        //this->tot_norms(i) = sqrt(ex_i.getSquareNorm());
    //}

    //double time = timer.getWallTime();
    //TelePrompter::printTree(0, "Hartree-Fock exchange", maxNodes, time);
}

int ExchangePotential::calcInternal(int i) {
    NOT_IMPLEMENTED_ABORT;
    //int maxNodes = 0;
    //Orbital &phi_i = this->orbitals_0->getOrbital(i);

    //double prec = getScaledPrecision(i, i);
    //prec = min(prec, 1.0e-1);
    //println(2, "\n [" << i << "," << i << "]");

    //Orbital *phi_ii = new Orbital(phi_i);
    //{ // compute phi_ii = phi_i^dag * phi_i
        //Timer timer;
        //this->mult.setPrecision(prec);
        //this->mult.adjoint(*phi_ii, 1.0, phi_i, phi_i);
        //int nNodes = phi_ii->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Multiply", nNodes, time);
    //}

    //Orbital *V_ii = new Orbital(phi_i);
    //{ // compute V_ii = P[phi_ii]
        //Timer timer;
        //this->poisson.setPrecision(prec);
        //if (phi_ii->hasReal()) {
            //V_ii->real = this->grid();
            //this->poisson(*V_ii->real, *phi_ii->real);
        //}
        //if (phi_ii->hasImag()) {
            //V_ii->imag = this->grid();
            //this->poisson(*V_ii->imag, *phi_ii->imag);
        //}
        //int nNodes = V_ii->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Applying Poisson", nNodes, time);
    //}
    //if (phi_ii != 0) delete phi_ii;

    //Orbital *phi_iii = new Orbital(phi_i);
    //{ // compute phi_iii = phi_i * V_ii
        //Timer timer;
        //double fac = -(this->x_factor/phi_i.getSquareNorm());
        //this->mult(*phi_iii, fac, phi_i, *V_ii);
        //this->part_norms(i,i) = sqrt(phi_iii->getSquareNorm());
        //int nNodes = phi_iii->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Multiply", nNodes, time);
    //}
    //if (V_ii != 0) delete V_ii;

    //{ // compute x_i += phi_iii
        //Timer timer;
        //Orbital &oldEx_i = this->exchange_0.getOrbital(i);
        //Orbital *newEx_i = new Orbital(oldEx_i);
        //this->add(*newEx_i, 1.0, oldEx_i, 1.0, *phi_iii);
        //int nNodes = newEx_i->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Adding exchange", nNodes, time);
        //this->exchange_0.replaceOrbital(i, &newEx_i);
    //}
    //if (phi_iii != 0) delete phi_iii;
    //return maxNodes;
}

int ExchangePotential::calcInternal(int i, int j) {
    NOT_IMPLEMENTED_ABORT;
    //int maxNodes = 0;

    //Orbital &phi_i = this->orbitals_0->getOrbital(i);
    //Orbital &phi_j = this->orbitals_0->getOrbital(j);

    //double i_factor = phi_i.getExchangeFactor(phi_j);
    //double j_factor = phi_j.getExchangeFactor(phi_i);
    //if (i_factor < MachineZero or j_factor < MachineZero) {
        //this->part_norms(i,j) = 0.0;
        //return 0;
    //}

    //printout(2, "\n [" << i << "," << j << "]");
    //printout(2, "   [" << j << "," << i << "]" << endl);
    //double prec = getScaledPrecision(i, j);
    //if (prec > 1.0e00) {
        //return 0;
    //}

    //Orbital *phi_ij = new Orbital(phi_i);
    //{ // compute phi_ij = phi_i^dag * phi_j
        //Timer timer;
        //if (phi_i.hasImag()) NOT_IMPLEMENTED_ABORT;
        //this->mult.adjoint(*phi_ij, 1.0, phi_i, phi_j);
        //int nNodes = phi_ij->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Multiply", nNodes, time);
    //}

    //Orbital *V_ij = new Orbital(phi_i);
    //{ // compute V_ij = P[phi_ij]
        //Timer timer;
        //prec = min(prec, 1.0e-1);
        //this->poisson.setPrecision(prec);
        //if (phi_ij->hasReal()) {
            //V_ij->real = this->grid();
            //this->poisson(*V_ij->real, *phi_ij->real);
        //}
        //if (phi_ij->hasImag()) {
            //V_ij->imag = this->grid();
            //this->poisson(*V_ij->imag, *phi_ij->imag);
        //}
        //int nNodes = V_ij->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Applying Poisson", nNodes, time);
    //}
    //if (phi_ij != 0) delete phi_ij;

    //Orbital *phi_iij = new Orbital(phi_i);
    //{ // compute phi_iij = phi_i * V_ij
        //Timer timer;
        //double fac = -(this->x_factor/phi_i.getSquareNorm());
        //this->mult(*phi_iij, fac, phi_i, *V_ij);
        //this->part_norms(i,j) = sqrt(phi_iij->getSquareNorm());
        //int nNodes = phi_iij->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Multiply", nNodes, time);
    //}

    //Orbital *phi_jij = new Orbital(phi_i);
    //{ // compute phi_jij = phi_j * V_ij
        //Timer timer;
        //double fac = -(this->x_factor/phi_j.getSquareNorm());
        //this->mult(*phi_jij, fac, phi_j, *V_ij);
        //this->part_norms(j,i) = sqrt(phi_jij->getSquareNorm());
        //int nNodes = phi_jij->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Multiply", nNodes, time);
    //}
    //if (V_ij != 0) delete V_ij;

    //{ // compute x_i += phi_jij
        //Timer timer;
        //Orbital &oldEx_i = this->exchange_0.getOrbital(i);
        //Orbital *newEx_i = new Orbital(oldEx_i);
        //this->add(*newEx_i, 1.0, oldEx_i, i_factor, *phi_jij);
        //int nNodes = newEx_i->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Adding exchange", nNodes, time);
        //this->exchange_0.replaceOrbital(i, &newEx_i);
    //}
    //if (phi_jij != 0) delete phi_jij;

    //{ // compute x_j += phi_iij
        //Timer timer;
        //Orbital &oldEx_j = this->exchange_0.getOrbital(j);
        //Orbital *newEx_j = new Orbital(oldEx_j);
        //this->add(*newEx_j, 1.0, oldEx_j, j_factor, *phi_iij);
        //int nNodes = newEx_j->getNNodes();
        //double time = timer.getWallTime();
        //maxNodes = max(maxNodes, nNodes);
        //TelePrompter::printTree(2, "Adding exchange", nNodes, time);
        //this->exchange_0.replaceOrbital(j, &newEx_j);
    //}
    //if (phi_iij != 0) delete phi_iij;
    //return maxNodes;
}

Orbital* ExchangePotential::testPreComputed(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
//    for (int i = 0; i < this->orbitals_0->size(); i++) {
//        Orbital &phi_i  = this->orbitals_0->getOrbital(i);
//        Orbital *ex_i = this->exchange_0.getOrbitalPtr(i);
//        if (&phi_i == &phi_p and ex_i != 0) {
//            Orbital *result = new Orbital(phi_p);
//            *result = *ex_i;
//            return result;
//        }
//    }
//    return 0;
}
