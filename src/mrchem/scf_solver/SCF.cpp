#include "SCF.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

SCF::SCF(const MultiResolutionAnalysis<3> &mra, HelmholtzOperatorSet &h)
    : nIter(0),
      maxIter(-1),
      iterPerRotation(0),
      orbThrs(-1.0),
      propThrs(-1.0),
      helmholtz(&h),
      add(mra, -1.0),
      MRA(mra) {
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

    this->add.setPrecision(this->orbPrec[0]);

    TelePrompter::printSeparator(0, '=');
    TelePrompter::printDouble(0, "Current precision", this->orbPrec[0]);
    TelePrompter::printSeparator(0, '-');
    TelePrompter::printDouble(0, "Orbital threshold", this->orbThrs);
    TelePrompter::printDouble(0, "Property threshold", this->propThrs);
    TelePrompter::printSeparator(0, '=', 2);
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

void SCF::printOrbitals(const MatrixXd &F, const OrbitalVector &phi) const {
    TelePrompter::printHeader(0, "Orbitals");
    println(0, " Orb    F(i,i)        Error         nNodes  Spin  Occ  Done ");
    TelePrompter::printSeparator(0, '-');
    int oldprec = TelePrompter::setPrecision(5);

    for (int i = 0; i < phi.size(); i++) {
        const Orbital &phi_i = phi.getOrbital(i);
        printout(0, setw(3) << i);
        printout(0, " " << setw(13) << F(i,i));
        printout(0, " " << setw(13) << phi_i.getError());
        printout(0, " " << setw(10) << phi_i.getNNodes());
        printout(0, setw(5) << phi_i.printSpin());
        printout(0, setw(5) << phi_i.getOccupancy());
        printout(0, setw(5) << phi_i.isConverged(getOrbitalThreshold()) << endl);
    }
    TelePrompter::printSeparator(0, '-');
    printout(0, " Total error:                    ");
    printout(0, setw(19) << phi.calcTotalError() << "  ");
    printout(0, setw(3) << phi.isConverged(getOrbitalThreshold()) << endl);
    TelePrompter::printSeparator(0, '=', 2);
    TelePrompter::setPrecision(oldprec);
}

void SCF::printConvergence(bool converged) const {
    int iter = this->orbError.size();
    int oldPrec = TelePrompter::getPrecision();
    TelePrompter::printHeader(0, "Convergence rate");
    println(0,"Iter   Orb Error       Property                     Update  ");
    TelePrompter::printSeparator(0, '-');
    for (int i = 0; i < iter; i++) {
        double prop_i = this->property[i];
        double propDiff = getUpdate(this->property, i+1, true);
        printout(0, setw(3) << i+1);
        TelePrompter::setPrecision(5);
        printout(0, setw(15) << this->orbError[i]);
        TelePrompter::setPrecision(15);
        printout(0, setw(26) << prop_i);
        TelePrompter::setPrecision(5);
        printout(0, setw(16) << propDiff);
        printout(0, endl);
    }
    TelePrompter::setPrecision(oldPrec);
    TelePrompter::printSeparator(0, '-');
    if (converged) {
        println(0,"                      SCF converged!!!                      ");
    } else {
        println(0,"                   SCF did NOT converge!!!                  ");
    }
    TelePrompter::printSeparator(0, '=', 2);
}

void SCF::printCycle() const {
    printout(0, endl << endl);
    printout(0, "#######################");
    printout(0, " SCF cycle " << setw(2) << this->nIter << " ");
    printout(0, "#######################");
    printout(0, endl << endl << endl);
}

void SCF::printTimer(double t) const {
    int oldPrec = TelePrompter::setPrecision(5);
    printout(0, endl << endl);
    printout(0, "################");
    printout(0, " Elapsed time:  " << t << " ");
    printout(0, "################");
    printout(0, endl << endl << endl);
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

/** Computes a new set of orbitals by application of the Helmholtz operator
 *
 * Requires orbitals in the orbitals pointer as well as corresponding operators
 * in the Helmholtz set and produces new non-orthogonal orbitals which are
 * stored locally in the newOrbs pointer. Does not compute orbital differences,
 * so the orbUpdates set is empty after this routine. Errors are estimated by
 * the norms of the new orbitals (deviation from one).
 */
void SCF::applyHelmholtzOperators(OrbitalVector &phi_np1,
                                  MatrixXd &F_n,
                                  OrbitalVector &phi_n,
                                  bool adjoint) {
    Timer timer;
    timer.restart();

    TelePrompter::printHeader(0, "Applying Helmholtz Operators");
    println(0, " Orb    OrbNorm       NormDiff       nNodes         Timing   ");
    TelePrompter::printSeparator(0, '-');
    int oldprec = TelePrompter::setPrecision(5);

    HelmholtzOperatorSet &H = *this->helmholtz;
    H.setPrecision(getOrbitalPrecision());

    phi_np1.clear();
    for (int i = 0; i < phi_n.size(); i++) {
        Timer timer;
        timer.restart();
        Orbital &nPhi_i = phi_n.getOrbital(i);
        Orbital &np1Phi_i = phi_np1.getOrbital(i);

        Orbital *arg_i = getHelmholtzArgument(i, F_n, phi_n, adjoint);
        H(i, np1Phi_i, *arg_i);
        delete arg_i;

        int nNodes = np1Phi_i.getNNodes();
        double norm_n = sqrt(nPhi_i.getSquareNorm());
        double norm_np1 = sqrt(np1Phi_i.getSquareNorm());
        double dNorm_n = fabs(norm_np1-norm_n);

        printout(0, setw(3) << i);
        printout(0, " " << setw(13) << norm_np1);
        printout(0, " " << setw(13) << dNorm_n);
        printout(0, " " << setw(9) << nNodes);
        printout(0, setw(18) << timer.getWallTime() << endl);
    }
    TelePrompter::printFooter(0, timer, 2);
    TelePrompter::setPrecision(oldprec);
}

Orbital* SCF::calcMatrixPart(int i, MatrixXd &M, OrbitalVector &phi) {
    vector<double> coefs;
    vector<Orbital *> orbs;

    int nOrbs = phi.size();
    for (int j = 0; j < nOrbs; j++) {
        double coef = M(i,j);
        // Linear scaling screening inserted here
        if (fabs(coef) > MachineZero) {
            Orbital &phi_j = phi.getOrbital(j);
            double norm_j = sqrt(phi_j.getSquareNorm());
            if (norm_j > 0.01*getOrbitalPrecision()) {
                coefs.push_back(coef);
                orbs.push_back(&phi_j);
            }
        }
    }

    Orbital &phi_i = phi.getOrbital(i);
    Orbital *result = new Orbital(phi_i);
    if (orbs.size() > 0) {
        Timer timer;
        timer.restart();
        this->add(*result, coefs, orbs);
        double time = timer.getWallTime();
        int nNodes = result->getNNodes();
        TelePrompter::printTree(2, "Added matrix part", nNodes, time);
    }
    return result;
}
