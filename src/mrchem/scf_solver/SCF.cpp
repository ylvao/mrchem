#include "SCF.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

SCF::SCF(HelmholtzOperatorSet &h)
    : nIter(0),
      maxIter(-1),
      iterPerRotation(0),
      orbThrs(-1.0),
      propThrs(-1.0),
      helmholtz(&h) {
    this->orbPrec[0] = -1.0;
    this->orbPrec[1] = -1.0;
    this->orbPrec[2] = -1.0;
}

SCF::~SCF() {
    this->helmholtz = 0;
}

void SCF::setThreshold(double orb_thrs, double prop_thrs) {
    this->orbThrs = orb_thrs;
    this->propThrs = prop_thrs;
}

void SCF::setOrbitalPrec(double init, double final) {
    this->orbPrec[0] = init;
    this->orbPrec[1] = init;
    this->orbPrec[2] = final;
}

void SCF::adjustPrecision(double error) {
    if (this->orbPrec[0] > 0.0 ) this->orbPrec[0] *= 0.75;
    this->orbPrec[0] = min(10.0*error*error, this->orbPrec[0]);
    this->orbPrec[0] = max(this->orbPrec[0], this->orbPrec[2]);

    printout(0, endl);
    TelePrompter::printDouble(0, "Current precision", this->orbPrec[0]);
    printout(0, endl);
    TelePrompter::printDouble(0, "Orbital threshold", this->orbThrs);
    TelePrompter::printDouble(0, "Property threshold", this->propThrs);
    printout(0, endl);
}

void SCF::resetPrecision() {
    this->orbPrec[0] = this->orbPrec[1];
}

bool SCF::needLocalization() const {
    if (this->iterPerRotation <= 0) {
        return false;
    }
    if (this->nIter <= 2) {
        return true;
    }
    if (this->nIter%this->iterPerRotation == 0) {
        return true;
    }
    return false;
}

bool SCF::needDiagonalization() const {
    if (this->iterPerRotation >= 0) {
        return false;
    }
    if (this->nIter <= 2) {
        return true;
    }
    if (this->nIter%this->iterPerRotation == 0) {
        return true;
    }
    return false;
}

void SCF::printUpdate(const string &name, double P, double dP) const {
    int oldPrec = TelePrompter::setPrecision(15);
    double p = 1.0;
    if (fabs(P) > MachineZero) {
        p = P;
    }
    bool done = fabs(dP/p) < getPropertyThreshold();
    printout(0, name);
    printout(0, setw(25) << P);
    TelePrompter::setPrecision(5);
    printout(0, setw(17) << dP);
    println(0, setw(5) << done);
    TelePrompter::setPrecision(oldPrec);
}

double SCF::getUpdate(const vector<double> &vec, int i, bool absPrec) const {
    if (i < 1 or i > vec.size()) MSG_ERROR("Invalid argument");
    double E_i = vec[i-1];
    double E_im1 = 0.0;
    if (i > 1) {
        E_im1 = vec[i-2];
    }
    double E_diff = E_i - E_im1;
    if (not absPrec and fabs(E_i) > MachineZero) {
        E_diff *= 1.0/E_i;
    }
    return E_diff;
}

void SCF::printOrbitals(const MatrixXd &f_mat, const OrbitalVector &phi) const {
    NOT_IMPLEMENTED_ABORT;
//    println(0, endl);
//    println(0, "------------------------------------------------------------");
//    println(0, "Orb    F(i,i)         Error         nNodes  Spin  Occ  Done ");
//    println(0, "------------------------------------------------------------");

//    int nOrbs = phi.size();
//    for (int i = 0; i < nOrbs; i++) {
//        const Orbital &phi_i = phi.getOrbital(i);

//        char sp = 'u';
//        if (phi_i.getSpin() == Orbital::Alpha) sp = 'a';
//        if (phi_i.getSpin() == Orbital::Beta) sp = 'b';

//        printout(0, setw(3) << i);
//        printout(0, " " << setw(13) << f_mat(i,i));
//        printout(0, " " << setw(13) << phi_i.getError());
//        printout(0, " " << setw(10) << phi_i.getNNodes());
//        printout(0, setw(5) << sp);
//        printout(0, setw(5) << phi_i.getOccupancy());
//        printout(0, setw(5) << phi_i.isConverged(getOrbitalThreshold()) << endl);
//    }
//    printout(0, "------------------------------");
//    println(0,  "------------------------------");
//    printout(0, "Total error:                     ");
//    printout(0, setw(16) << phi.calcTotalError() << "     ");
//    printout(0, setw(3) << phi.isConverged(getOrbitalThreshold()) << endl);
//    printout(0, "------------------------------");
//    println(0,  "------------------------------");
//    printout(0, "\n");
}

void SCF::printOrbitals(const OrbitalVector &phi) const {
    NOT_IMPLEMENTED_ABORT;
//    TelePrompter::printHeader(0, "Orbital Vector");
//    println(0, "Orb    Norm           Error         nNodes  Spin  Occ  Done ");
//    TelePrompter::printSeparator(0, '-');

//    int nOrbs = phi.size();
//    for (int i = 0; i < nOrbs; i++) {
//        const Orbital &phi_i = phi.getOrbital(i);

//        char sp = 'u';
//        if (phi_i.getSpin() == Orbital::Alpha) sp = 'a';
//        if (phi_i.getSpin() == Orbital::Beta) sp = 'b';

//        printout(0, setw(3) << i);
//        printout(0, " " << setw(13) << sqrt(phi_i.getSquareNorm()));
//        printout(0, " " << setw(13) << phi_i.getError());
//        printout(0, " " << setw(10) << phi_i.getNNodes());
//        printout(0, setw(5) << sp);
//        printout(0, setw(5) << phi_i.getOccupancy());
//        printout(0, setw(5) << phi_i.isConverged(getOrbitalThreshold()) << endl);
//    }
//    TelePrompter::printSeparator(0, "-");
//    printout(0, "Total error:                     ");
//    printout(0, setw(16) << phi.calcTotalError() << "     ");
//    printout(0, setw(3) << phi.isConverged(getOrbitalThreshold()) << endl);
//    TelePrompter::printSeparator(0, "-", 2);
}

void SCF::printConvergence(bool converged) const {
    NOT_IMPLEMENTED_ABORT;
//    int iter = this->orbError.size();
//    int oldPrec = TelePrompter::getPrecision();
//    println(0,"                                                            ");
//    println(0,"                                                            ");
//    println(0,"==================== Convergence rate ======================");
//    println(0,"                                                            ");
//    if (converged) {
//        println(0, " SCF converged in " << iter << " iterations!");
//    } else {
//        println(0, " SCF did NOT converge in " << iter << " iterations!");
//    }
//    println(0,"                                                            ");
//    println(0,"------------------------------------------------------------");
//    println(0,"Iter   Orb Error       Property                     Update  ");
//    println(0,"------------------------------------------------------------");
//    for (int i = 0; i < iter; i++) {
//        double prop_i = this->property[i];
//        double propDiff = getUpdate(this->property, i+1, true);
//        printout(0, setw(3) << i+1);
//        TelePrompter::setPrecision(5);
//        printout(0, setw(15) << this->orbError[i]);
//        TelePrompter::setPrecision(15);
//        printout(0, setw(26) << prop_i);
//        TelePrompter::setPrecision(5);
//        printout(0, setw(16) << propDiff);
//        printout(0, endl);
//    }
//    TelePrompter::setPrecision(oldPrec);
//    println(0,"============================================================");
//    println(0,"                                                            ");
//    println(0,"                                                            ");
}

void SCF::printCycle() const {
    printout(0, endl << endl);
    printout(0, "#######################");
    printout(0, " SCF cycle " << setw(2) << this->nIter << " ");
    printout(0, "#######################");
    printout(0, endl << endl);
}

void SCF::printTimer(double t) const {
    int oldPrec = TelePrompter::setPrecision(5);
    printout(0, endl << endl);
    printout(0, "################");
    printout(0, " Elapsed time:  " << t << " ");
    printout(0, "################");
    printout(0, endl << endl);
    TelePrompter::setPrecision(oldPrec);
}

void SCF::printMatrix(int level, const MatrixXd &M, const char &name, int pr) const {
    int oldPrec = TelePrompter::setPrecision(pr);
    printout(level, endl);
    printout(level, "----------------------------- ");
    printout(level, name);
    printout(level, " ----------------------------");
    printout(level, endl);
    printout(level, M);
    printout(level, endl);
    printout(level, "------------------------------");
    printout(level, "------------------------------");
    printout(level, endl);
    printout(level, endl);
    TelePrompter::setPrecision(oldPrec);
}

bool SCF::accelerate(Accelerator *acc, OrbitalVector *phi, OrbitalVector *d_phi,
                     MatrixXd *f_mat, MatrixXd *df_mat) {
    NOT_IMPLEMENTED_ABORT;
//    if (acc != 0) {
//        if (phi == 0) MSG_ERROR("Uninitialized orbitals");
//        if (d_phi == 0) MSG_ERROR("Uninitialized orbitals");
//        acc->pushBack(*phi, *d_phi, f_mat, df_mat);
//        acc->calcUpdates(*phi, *d_phi, f_mat, df_mat);
//        return true;
//    }
//    return false;
}

/** Computes a new set of orbitals by application of the Helmholtz operator
 *
 * Requires orbitals in the orbitals pointer as well as corresponding operators
 * in the Helmholtz set and produces new non-orthogonal orbitals which are
 * stored locally in the newOrbs pointer. Does not compute orbital differences,
 * so the orbUpdates set is empty after this routine. Errors are estimated by
 * the norms of the new orbitals (deviation from one).
 */
void SCF::applyHelmholtzOperators(OrbitalVector &phi_np1,
                                  Eigen::MatrixXd &f_mat_n,
                                  OrbitalVector &phi_n,
                                  bool adjoint) {
    NOT_IMPLEMENTED_ABORT;
//    Timer timer;
//    timer.restart();

//    TelePrompter::printHeader(0, "Calculating Orbital Updates");
//    println(0, "Orb    OrbNorm       NormDiff       nNodes         Timing   ");

//    phi_np1.clear();
//    for (int i = 0; i < phi_n.size(); i++) {
//        Timer timer;
//        timer.restart();
//        Orbital &nPhi_i = phi_n.getOrbital(i);
//        double norm_n = nPhi_i.getSquareNorm();
//        if (norm_n < 0.0) {
//            norm_n = 0.0;
//        } else {
//            norm_n = sqrt(norm_n);
//        }

//        HelmholtzOperator &H = this->helmholtz->getOperator(i);
//        Orbital *greenArgument = getHelmholtzArgument(i, phi_n, f_mat_n, adjoint);

//        Orbital *np1Phi_i = new Orbital(nPhi_i);
//        np1Phi_i->setRelPrec(getOrbitalPrecision());
//        np1Phi_i->setAbsPrec(true);
//        greenOper.apply(*np1Phi_i, *greenArgument);
//        delete greenArgument;

//        int nNodes = np1Phi_i->getNNodes();
//        double norm_np1 = sqrt(np1Phi_i->getSquareNorm());
//        double dNorm_n = fabs(norm_np1-norm_n);

//        printout(0, setw(3) << i);
//        printout(0, " " << setw(13) << norm_np1);
//        printout(0, " " << setw(13) << dNorm_n);
//        printout(0, " " << setw(10) << nNodes);
//        printout(0, setw(18) << timer.elapsed() << endl);

//        phi_np1.replaceOrbital(i, &np1Phi_i);
//    }
//    TelePrompter::printFooter(0, timer, 2);
}

Orbital* SCF::calcMatrixPart(int i, MatrixXd &M, OrbitalVector &phi) {
    NOT_IMPLEMENTED_ABORT;
//    boost::timer timer;
//    vector<double> expCoefs;
//    vector<FunctionTree<3> *> expOrbs;

//    int nOrbs = phi.size();
//    for (int j = 0; j < nOrbs; j++) {
//        double coef = M(i,j);
//        // Linear scaling screening inserted here
//        if (fabs(coef) > MachineZero) {
//            Orbital &orb_j = phi.getOrbital(j);
//            double norm_j = orb_j.getSquareNorm();
//            if (norm_j > 0.0) {
//                norm_j = sqrt(norm_j);
//            } else {
//                norm_j = 0.0;
//            }
//            if (norm_j > 0.01*getOrbitalPrecision()) {
//                expCoefs.push_back(coef);
//                expOrbs.push_back(&orb_j);
//            }
//        }
//    }

//    Orbital *orb = 0;
//    if (expOrbs.size() > 0) {
//        timer.restart();
//        Orbital &orb_i = phi.getOrbital(i);
//        orb = new Orbital(orb_i);
//        orb->setRelPrec(getOrbitalPrecision());
//        orb->setAbsPrec(true);
//        orb->add(expCoefs, expOrbs, -1);
//        double time = timer.elapsed();
//        int nNodes = orb->getNNodes();
//        TelePrompter::printTree(2, "Added matrix part", nNodes, time);
//    }
//    return orb;
}

