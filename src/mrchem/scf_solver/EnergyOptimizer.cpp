#include "EnergyOptimizer.h"

using namespace std;
using namespace Eigen;

EnergyOptimizer::EnergyOptimizer(const MultiResolutionAnalysis<3> &mra,
                                 HelmholtzOperatorSet &h)
        : GroundStateSolver(mra, h),
          fOper_np1(0) {
}

void EnergyOptimizer::setup(FockOperator &fock, OrbitalVector &phi, MatrixXd &F,
               FockOperator &fock_np1, OrbitalVector &phi_np1) {
    NOT_IMPLEMENTED_ABORT;
//    CoulombOperator *j_np1 = 0;
//    ExchangeOperator *k_np1 = 0;
//    XCOperator *xc_np1 = 0;

//    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
//    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
//    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
//    XCOperator *xc_n = this->fOper_n->getXCOperator();

//    double prec = this->orbPrec[2];
//    if (j_n != 0) j_np1 = new CoulombPotential(prec, this->MRA, *this->orbitals_np1);
//    if (k_n != 0) k_np1 = new ExchangePotential(prec, this->MRA, *this->orbitals_np1);
//    if (xc_n != 0) xc_np1 = new XCPotential(prec, this->MRA, *this->orbitals_np1);

//    this->fOper_np1 = new FockOperator(this->MRA, 0, v_n, j_np1, k_np1, xc_np1);
}

void EnergyOptimizer::clear() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->fOper_np1 != 0) {
//        CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
//        ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
//        XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

//        if (j_np1 != 0) delete j_np1;
//        if (k_np1 != 0) delete k_np1;
//        if (xc_np1 != 0) delete xc_np1;

//        delete this->fOper_np1;
//    }
//    this->fOper_np1 = 0;
}

bool EnergyOptimizer::optimize() {
    NOT_IMPLEMENTED_ABORT;
//    boost::timer scfTimer;
//    double err_p = 1.0;
//    double err_o = this->phi_n->getErrors().maxCoeff();
//    double err_t = 1.0;

//    bool first = true;
//    bool converged = false;
//    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
//        scfTimer.restart();
//        printCycle();
//        adjustPrecision(err_o);

//        printMatrix(1, this->phi_n->calcOverlapMatrix(), 'S');
//        printMatrix(1, *this->fMat_n, 'F');

//        this->fOper_n->setup(getOrbitalPrecision());
//        double prop = calcProperty();
//        this->property.push_back(prop);

//        this->helmholtz->initialize(this->fMat_n->diagonal());
//        applyHelmholtzOperators(*this->phi_np1, *this->fMat_n, *this->phi_n);
//        calcOrbitalUpdates();
//        calcFockMatrixUpdate();
//        printTreeSizes();
//        this->dPhi_n->clear();

//        *this->fMat_np1 = *this->fMat_n + *this->dfMat_n;

//        printMatrix(1, *this->fMat_np1, 'F');
//        rotate(*this->phi_np1, *this->fMat_np1);
//        printMatrix(1, *this->fMat_np1, 'F');

//        calcOrbitalUpdates();
//        *this->dfMat_n = *this->fMat_np1 - *this->fMat_n;
//        this->phi_np1->clear();

//        err_p = calcPropertyError();
//        err_o = calcOrbitalError();
//        err_t = calcTotalError();
//        this->orbError.push_back(err_t);

//        this->phi_n->addInPlace(*this->dPhi_n);
//        *this->fMat_n += *this->dfMat_n;

//        this->phi_n->orthonormalize(getOrbitalPrecision(), this->fMat_n);
//        this->phi_n->setErrors(this->dPhi_n->getNorms());

//        this->fOper_n->clear();
//        this->dPhi_n->clear();

//        printOrbitals(*this->fMat_n, *this->phi_n);
//        printProperty();
//        printTimer(scfTimer.elapsed());

//        if (err_p < getPropertyThreshold() and not first) {
//            converged = true;
//            break;
//        }
//        first = false;
//    }
//    return converged;
}

void EnergyOptimizer::calcFockMatrixUpdate() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");
//    boost::timer timer, totTimer;
//    println(0,"                                                            ");
//    println(0,"=============== Computing Fock Matrix update ===============");
//    println(0,"                                                            ");

//    timer.restart();
//    MatrixXd dS_1 = this->dPhi_n->calcOverlapMatrix(*this->phi_n);
//    MatrixXd dS_2 = this->phi_np1->calcOverlapMatrix(*this->dPhi_n);
//    double time_s = timer.elapsed();
//    TelePrompter::printDouble(0, "Overlap matrices", time_s);

//    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
//    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
//    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
//    XCOperator *xc_n = this->fOper_n->getXCOperator();

//    // Nuclear potential matrix is computed explicitly
//    timer.restart();
//    MatrixXd dV_n = (*v_n)(*this->phi_np1, *this->dPhi_n);
//    double time_v = timer.elapsed();
//    TelePrompter::printDouble(0, "Nuclear potential matrix", time_v);

//    timer.restart();
//    FockOperator f_n(0, 0, j_n, k_n, xc_n);
//    MatrixXd F_n = f_n(*this->phi_np1, *this->phi_n);
//    double time_f_n = timer.elapsed();
//    TelePrompter::printDouble(0, "Fock matrix n", time_f_n);
//    this->fOper_n->clear();

//    // The n+1 Fock operator needs orthonormalized orbitals
//    this->phi_np1->orthonormalize(getOrbitalPrecision());

//    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
//    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
//    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

//    println(0,"                                                            ");
//    // Do not setup exchange, it must be applied on the fly anyway
//    if (j_np1 != 0) j_np1->setup(getOrbitalPrecision());
//    if (k_np1 != 0) j_np1->QMOperator::setup(getOrbitalPrecision());
//    if (xc_np1 != 0) xc_np1->setup(getOrbitalPrecision());
//    println(0,"                                                            ");

//    // Computing potential matrix excluding nuclear part
//    timer.restart();
//    FockOperator f_np1(0, 0, j_np1, k_np1, xc_np1);
//    MatrixXd F_1 = f_np1(*this->phi_n, *this->phi_n);
//    MatrixXd F_2 = f_np1(*this->phi_n, *this->dPhi_n);
//    MatrixXd F_np1 = F_1 + F_2 + F_2.transpose();
//    //MatrixXd F_3 = f_np1(*this->dPhi_n, *this->phi_n);
//    //MatrixXd F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
//    //MatrixXd F_np1 = F_1 + F_2 + F_3 + F_4;
//    double time_f_np1 = timer.elapsed();
//    TelePrompter::printDouble(0, "Fock matrix n+1", time_f_np1);

//    f_np1.clear();
//    this->phi_np1->clear();
//    this->phi_np1->add(1.0, *this->phi_n, 1.0, *this->dPhi_n);

//    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
//    MatrixXd dF_n = F_np1 - F_n;
//    MatrixXd dF_1 = dS_1*(*this->fMat_n);
//    MatrixXd dF_2 = dS_2*L;

//    // Adding up the pieces
//    *this->dfMat_n = dF_1 + dF_2 + dV_n + dF_n;

//    //Symmetrizing
//    MatrixXd sym = *this->dfMat_n + this->dfMat_n->transpose();
//    *this->dfMat_n = 0.5 * sym;

//    double t = totTimer.elapsed();
//    println(0, "                                                            ");
//    println(0, "================ Elapsed time: " << t << " =================");
//    println(0, "                                                            ");
//    println(0, "                                                            ");
}
