#include <Eigen/Eigenvalues>

#include "OrbitalVector.h"
#include "Orbital.h"
//#include "PositionFunction.h"
//#include "HydrogenicFunction.h"
//#include "NonlinearMaximizer.h"
//#include "eigen_disable_warnings.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

/** OrbitalVector constructor
 *
 * New orbitals are constructed with double
 * occupancy and undefined spin
 *
 * All orbital functions are uninitialized.
 */
OrbitalVector::OrbitalVector(int n_orbs) {
    push_back(n_orbs, 2, Paired);
}

/** OrbitalVector constructor
 *
 * New orbitals are constructed with single
 * occupancy orbitals and alpha/beta spin
 *
 * All orbital functions are uninitialized.
 */
OrbitalVector::OrbitalVector(int n_alpha, int n_beta) {
    push_back(n_alpha, 1, Alpha);
    push_back(n_beta, 1, Beta);
}

/** OrbitalVector constructor
 *
 * ne is number of electrons.
 * mult is spin mutiplicity.
 * rest is spin restricted.
 *
 * All orbital functions are uninitialized.
 */
OrbitalVector::OrbitalVector(int ne, int mult, bool rest) {
    int de = ne - (mult - 1);
    if (de%2 != 0)  MSG_ERROR("Invalid multiplicity");

    int na, nb, nd;
    if (rest) {
        na = mult - 1;
        nb = 0;
        nd = de/2;
    } else {
        na = de/2 + (mult - 1);
        nb = de/2;
        nd = 0;
    }

    push_back(nd, 2, Paired);
    push_back(na, 1, Alpha);
    push_back(nb, 1, Beta);
}

/** Copy constructor
 *
 * New orbitals are constructed with spin and occupancy
 * parameters (not function data) taken from the input set.
 *
 * All orbital functions are uninitialized.
 */
OrbitalVector::OrbitalVector(const OrbitalVector &orb_set) {
    for (int i = 0; i < orb_set.size(); i++) {
        const Orbital &orb_i = orb_set.getOrbital(i);
        Orbital *newOrb = new Orbital(orb_i);
        this->orbitals.push_back(newOrb);
    }
}

/** OrbitalVector destructor
 *
 * Deletes all orbitals in the vector
 */
OrbitalVector::~OrbitalVector() {
    for (int i = 0; i < this->size(); i++) {
        if (this->orbitals[i] != 0) {
            delete this->orbitals[i];
        }
    }
    this->orbitals.clear();
}

/** Clears each orbital in the vector
 *
 * Deletes the actual functions in the orbitals, keeps
 * the spin and occupancy.
 */
void OrbitalVector::clear() {
    for (int i = 0; i < this->size(); i++) {
        this->orbitals[i]->clear();
    }
}

/** Append orbital to this set
 *
 * n_orbs is number of new orbitals.
 * occ is occupancy of all new orbitals.
 * spin is the spin of all new orbitals.
 *
 * New orbitals are constructed with given spin and occupancy
 * parameters, as uninitialized functions.
 *
 * Any existing orbitals in the set are kept.
 */
void OrbitalVector::push_back(int n_orbs, int occ, int spin) {
    for (int i = 0; i < n_orbs; i++) {
        Orbital *orb = new Orbital(occ, spin);
        this->orbitals.push_back(orb);
    }
}

//void OrbitalVector::initialize(const Nuclei &nucs) {
//    NOT_IMPLEMENTED_ABORT;
//    this->clear();
//    this->orbitals.clear();
//    for (int i = 0; i < nucs.size(); i++) {
//        const Nucleus &nuc = *nucs[i];
//        int minOrbs = ceil(nuc.getElement().getZ()/2.0);
//        double Z = nuc.getCharge();
//        const double *R = nuc.getCoord();
//        int n = 1;
//        int nOrbs = 0;
//        while (nOrbs < minOrbs) {
//            for (int l = 0; l < n; l++) {
//                if (nOrbs >= minOrbs) continue;
//                int M = 2*l+1;
//                for (int m = 0; m < M; m++) {
//                    HydrogenicFunction h_func(n, l, m, Z, R);
//                    Orbital *h_orb = new Orbital(0, Orbital::Paired);
//                    h_orb->projectFunction(h_func);
//                    this->orbitals.push_back(h_orb);
//                    println(0, "Orb " << size() << " " << n << " " << l << " " << m << "   " << h_orb->getSquareNorm());
//                    nOrbs++;
//                }
//            }
//            n++;
//        }
//    }
//}

//OrbitalVector& OrbitalVector::operator=(const OrbitalVector &orb_set) {
//    NOT_IMPLEMENTED_ABORT;
//    if (this == &orb_set) {
//        return *this;
//    }
//    if (this->size() > orb_set.size()) MSG_ERROR("Size mismatch");

//    int nOrbs = this->size();
//    for (int i = 0; i < nOrbs; i++) {
//        const Orbital *thatOrb = orb_set.getOrbitalPtr(i);
//        Orbital *thisOrb = this->getOrbitalPtr(i);
//        if (thatOrb == 0) MSG_ERROR("Reading NULL orbital");
//        if (thisOrb == 0) MSG_ERROR("Orbital not initialized");
//        *thisOrb = *thatOrb;
//    }
//    return *this;
//}

//void OrbitalVector::writeOrbitals(const string &of) {
//    for (int i = 0; i < this->size(); i++) {
//        stringstream name;
//        name << "orbitals/" << of << "_";
//        if (i < 1000) name << "0";
//        if (i < 100) name << "0";
//        if (i < 10) name << "0";
//        name << i;
//        Orbital &orb = getOrbital(i);
//        orb.saveTree(name.str());
//    }
//}

//void OrbitalVector::readOrbitals(const string &of) {
//    for (int i = 0; i < this->size(); i++) {
//        stringstream name;
//        name << "orbitals/" << of << "_";
//        if (i < 1000) name << "0";
//        if (i < 100) name << "0";
//        if (i < 10) name << "0";
//        name << i;
//        Orbital &orb = getOrbital(i);
//        orb.loadTree(name.str());
//    }
//}

const Orbital& OrbitalVector::getOrbital(int i) const {
    if (this->orbitals[i] == 0) MSG_ERROR("Incomplete set");
    return *this->orbitals[i];
}

Orbital& OrbitalVector::getOrbital(int i) {
    if (this->orbitals[i] == 0) MSG_ERROR("Incomplete set");
    return *this->orbitals[i];
}

/** Returns the number of occupied orbitals */
int OrbitalVector::getNOccupied() const {
    int nOccupied = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getOccupancy() > 0) {
            nOccupied++;
        }
    }
    return nOccupied;
}

/** Returns the number of empty orbitals */
int OrbitalVector::getNEmpty() const {
    int nEmpty = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getOccupancy() == 0) {
            nEmpty++;
        }
    }
    return nEmpty;
}

/** Returns the number of singly occupied orbitals */
int OrbitalVector::getNSingly() const {
    int nSingly = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getOccupancy() == 1) {
            nSingly++;
        }
    }
    return nSingly;
}

/** Returns the number of paired orbitals */
int OrbitalVector::getNPaired() const {
    int nPaired = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getSpin() == Paired) {
            nPaired++;
        }
    }
    return nPaired;
}

/** Returns the number of alpha orbitals */
int OrbitalVector::getNAlpha() const {
    int nAlpha = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getSpin() == Alpha) {
            nAlpha++;
        }
    }
    return nAlpha;
}

/** Returns the number of beta orbitals */
int OrbitalVector::getNBeta() const {
    int nBeta = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getSpin() == Beta) {
            nBeta++;
        }
    }
    return nBeta;
}

/** Returns the number of doubly occupied orbitals */
int OrbitalVector::getNDoubly() const {
    int nDoubly = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (orb->getOccupancy() == 2) {
            nDoubly++;
        }
    }
    return nDoubly;
}

/** Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 */
int OrbitalVector::getNElectrons(int inpSpin) const {
    int nElectrons = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        int thisSpin = orb->getSpin();
        if (inpSpin == Paired) {
            nElectrons += orb->getOccupancy();
        } else if (inpSpin == Alpha) {
            if (thisSpin == Paired or thisSpin == Alpha) {
                nElectrons += 1;
            }
        } else if (inpSpin == Beta) {
            if (thisSpin == Paired or thisSpin == Beta) {
                nElectrons += 1;
            }
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

int OrbitalVector::getMultiplicity() const {
    NOT_IMPLEMENTED_ABORT;
}

bool OrbitalVector::isConverged(double prec) const {
    bool converged = true;
    for (int i = 0; i < this->size(); i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        if (not orb->isConverged(prec)) {
            converged = false;
        }
    }
    return converged;
}

double OrbitalVector::calcTotalError() const {
    const VectorXd &error = getErrors();
    return sqrt(error.dot(error));
}

/** Returns a vector containing the orbital errors */
VectorXd OrbitalVector::getErrors() const {
    int nOrbs = this->size();
    VectorXd errors = VectorXd::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        errors(i) = orb->getError();
    }
    return errors;
}

/** Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 */
void OrbitalVector::setErrors(const VectorXd &errors) {
    int nOrbs = this->size();
    if (nOrbs != errors.size()) {
        MSG_ERROR("Size mismatch");
    }
    for (int i = 0; i < nOrbs; i++) {
        Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        orb->setError(errors(i));
    }
}

/** Returns a vector containing the orbital spins
 */
VectorXi OrbitalVector::getSpins() const {
    int nOrbs = this->size();
    VectorXi spins = VectorXi::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        spins(i) = orb->getSpin();
    }
    return spins;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 */
void OrbitalVector::setSpins(const VectorXi &spins) {
    int nOrbs = this->size();
    if (nOrbs != spins.size()) {
        MSG_ERROR("Size mismatch");
    }
    for (int i = 0; i < nOrbs; i++) {
        Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        orb->setSpin(spins(i));
    }
}

/** Returns a vector containing the orbital occupancies
 */
VectorXi OrbitalVector::getOccupancies() const {
    int nOrbs = this->size();
    VectorXi occ = VectorXi::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        occ(i) = orb->getOccupancy();
    }
    return occ;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 */
void OrbitalVector::setOccupancies(const VectorXi &occ) {
    int nOrbs = this->size();
    if (nOrbs != occ.size()) {
        MSG_ERROR("Size mismatch");
    }
    for (int i = 0; i < nOrbs; i++) {
        Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        orb->setOccupancy(occ(i));
    }
}

/** Returns a vector containing the orbital square norms
 */
VectorXd OrbitalVector::getSquareNorms() const {
    int nOrbs = this->size();
    VectorXd norms = VectorXd::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        norms(i) = orb->getSquareNorm();
    }
    return norms;
}

/** Returns a vector containing the orbital norms
 */
VectorXd OrbitalVector::getNorms() const {
    int nOrbs = this->size();
    VectorXd norms = VectorXd::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital *orb = getOrbitalPtr(i);
        if (orb == 0) {
            continue;
        }
        norms(i) = sqrt(orb->getSquareNorm());
    }
    return norms;
}

//void OrbitalVector::replaceOrbitals(OrbitalVector &new_orbs) {
//    int nOrbs = this->size();
//    if (nOrbs != new_orbs.size()) {
//        MSG_ERROR("Size mismatch");
//    }
//    for (int i = 0; i < nOrbs; i++) {
//        Orbital *newOrb = new_orbs.getOrbitalPtr(i);
//        replaceOrbital(i, &newOrb);
//        new_orbs.orbitals[i] = 0;
//    }
//}

void OrbitalVector::replaceOrbital(int i, Orbital **orb) {
    if (i < 0 or i >= this->size()) {
        MSG_ERROR("Orbital index out of bounds");
    }
    if (this->orbitals[i] != 0) {
        delete this->orbitals[i];
    }
    this->orbitals[i] = *orb;
    *orb = 0;
}

//void OrbitalVector::rotate(double prec, const MatrixXd &U) {
//    boost::timer timer, totTimer;
//    totTimer.restart();

//    int nOrbs = this->size();
//    if (nOrbs < 1) {
//        MSG_ERROR("Cannot rotate empty set");
//        return;
//    }
//    if (nOrbs == 1) {
//        Orbital *orb = this->orbitals[0];
//        if (orb != 0) *orb *= U(0,0);
//        return;
//    }

//    double threshold;
//    if (prec > 0.0) {
//        threshold = 0.01*prec;
//    } else {
//        threshold = MachineZero;
//    }

//    VectorXd errors = this->getErrors();
//    OrbitalVector *newSet = new OrbitalVector("tmpSet", *this);

//    double time;
//    int nNodes;
//    int nSkip = 0;
//    int nDone = 0;
//    for (int i = 0; i < nOrbs; i++) {
//        vector<double> coefs;
//        vector<FunctionTree<3> *> funcs;

//        int n = 0;
//        Orbital *newOrb = 0;
//        for (int j = 0; j < nOrbs; j++) {
//            double coef = U(i,j);
//            Orbital *jOrb = getOrbitalPtr(j);

//            if (fabs(coef) <= threshold or jOrb == 0) {
//                nSkip++;
//                continue;
//            }

//            coefs.push_back(coef);
//            funcs.push_back(jOrb);

//            if (n == 0) {
//                //Copy parameters from first contributing orbital
//                newOrb = new Orbital(*jOrb);
//                if (prec > 0.0) newOrb->setRelPrec(prec);
//            } else {
//                if (newOrb->compareSpin(*jOrb) < 0) {
//                    MSG_WARN("Mixing incompatible functions with coef " << coef);
//                }
//            }
//            nDone++;
//            n++;
//        }
//        if (n > 0) {
//            timer.restart();
//            if (n < 10) {
//                newOrb->add(coefs, funcs, 0);
//                newOrb->crop(prec);
//            } else {
//                newOrb->add(coefs, funcs, -1);
//            }
//            time = timer.elapsed();
//            nNodes = newOrb->getNNodes();
//            TelePrompter::printTree(1, "Rotated orbital", nNodes, time);
//        }
//        newSet->replaceOrbital(i, &newOrb);
//    }
//    //println(0, "Rotated 	      " << setw(30) << nDone << " entries");
//    //println(0, "Skipped	      " << setw(30) << nSkip << " entries");
//    replaceOrbitals(*newSet);
//    delete newSet;

//    setErrors(errors);

//    printout(0, "Rotating orbitals                                ");
//    println(0, totTimer.elapsed());
//}

//void OrbitalVector::add(double a, OrbitalVector &set_a, double b, OrbitalVector &set_b) {
//    boost::timer timer;
//    double time;
//    int nNodes;
//    int nOrbs = this->size();
//    if (set_a.size() != nOrbs or set_b.size() != nOrbs) {
//        MSG_ERROR("Size mismatch");
//    }
//    this->clear();
//    for (int i = 0; i < nOrbs; i++) {
//        timer.restart();
//        Orbital &aOrb = set_a.getOrbital(i);
//        Orbital &bOrb = set_b.getOrbital(i);
//        aOrb.compare(bOrb);
//        Orbital *newOrb = new Orbital(aOrb);
//        newOrb->add(a, aOrb, b, bOrb);
//        time = timer.elapsed();
//        nNodes = newOrb->getNNodes();
//        TelePrompter::printTree(1, "Added orbital", nNodes, time);
//        replaceOrbital(i, &newOrb);
//    }
//}

//void OrbitalVector::addInPlace(OrbitalVector &func_set, double c) {
//    boost::timer timer;
//    double time;
//    int nNodes;
//    int nOrbs = this->size();
//    if (func_set.size() != nOrbs) {
//        MSG_ERROR("Size mismatch");
//    }

//    VectorXd errors = this->getErrors();

//    for (int i = 0; i < nOrbs; i++) {
//        timer.restart();
//        Orbital &oldOrb = this->getOrbital(i);
//        Orbital &diffOrb = func_set.getOrbital(i);
//        oldOrb.compare(diffOrb);
//        Orbital *newOrb = new Orbital(oldOrb);
//        newOrb->add(1.0, oldOrb, c, diffOrb);
//        time = timer.elapsed();
//        nNodes = newOrb->getNNodes();
//        TelePrompter::printTree(1, "Added orbital", nNodes, time);

//        this->replaceOrbital(i, &newOrb);
//    }

//    setErrors(errors);
//}

//void OrbitalVector::crop(double prec) {
//    boost::timer timer;
//    timer.restart();
//    println(0, "\n--------------------- Cropping orbitals --------------------");
//    int nOrbs = this->size();
//    for (int i = 0; i < nOrbs; i++) {
//        Orbital *orb = getOrbitalPtr(i);
//        if (orb == 0) continue;
//        int oldNodes = orb->getNNodes();
//        orb->crop(prec);
//        int newNodes = orb->getNNodes();
//        printout(0, setw(3) << i);
//        printout(0, setw(45) << oldNodes << "   -> ");
//        printout(0, setw(5) << newNodes);
//        printout(0, endl);
//    }
//    printout(0, "---------------- Elapsed time: " << timer.elapsed());
//    println(0, " -----------------\n");
//}

//MatrixXd OrbitalVector::getSpinMatrix() const {
//    int N = this->size();
//    MatrixXd S = MatrixXd::Zero(N,N);
//    int i = 0;
//    for (int n = 0; n < N; n++) {
//        const Orbital &phi_n = getOrbital(n);
//        if (phi_n.getSpin() == Orbital::Paired) {
//            S(i,n) = 1.0;
//            i++;
//        }
//    }
//    for (int n = 0; n < N; n++) {
//        const Orbital &phi_n = getOrbital(n);
//        if (phi_n.getSpin() == Orbital::Alpha) {
//            S(i,n) = 1.0;
//            i++;
//        }
//    }
//    for (int n = 0; n < N; n++) {
//        const Orbital &phi_n = getOrbital(n);
//        if (phi_n.getSpin() == Orbital::Beta) {
//            S(i,n) = 1.0;
//            i++;
//        }
//    }
//    return S;
//}

//void OrbitalVector::spinCleanMatrix(MatrixXd &M) const {
//    int nOrbs = this->size();
//    if (nOrbs != M.cols() or nOrbs != M.rows()) MSG_ERROR("Size mismatch");
//    for (int i = 0; i < nOrbs; i++) {
//        int spin_i = getOrbital(i).getSpin();
//        for (int j = 0; j < nOrbs; j++) {
//            int spin_j = getOrbital(j).getSpin();
//            if (spin_i == Orbital::Alpha and spin_j == Orbital::Beta) {
//                if (M(i,j) > MachineZero) MSG_WARN("Removing spin mixing");
//                M(i,j) = 0.0;
//            }
//            if (spin_i == Orbital::Beta and spin_j == Orbital::Alpha) {
//                if (M(i,j) > MachineZero) MSG_WARN("Removing spin mixing");
//                M(i,j) = 0.0;
//            }
//        }
//    }
//}

//void OrbitalVector::separateSpinMatrix(MatrixXd &T, MatrixXd &P, MatrixXd &A, MatrixXd &B) const {
//    MatrixXd S = getSpinMatrix();
//    T = S*T;

//    int N_p = getNPaired();
//    int N_a = getNAlpha();
//    int N_b = getNBeta();
//    P = T.block(0,0,N_p,N_p);
//    A = T.block(N_p,N_p,N_a,N_a);
//    B = T.block(N_p+N_a,N_p+N_a,N_b,N_b);
//}

//void OrbitalVector::collectSpinMatrix(MatrixXd &T, MatrixXd &P, MatrixXd &A, MatrixXd &B) const {
//    int rows = P.rows() + A.rows() + B.rows();
//    int cols = P.cols() + A.cols() + B.cols();
//    T = MatrixXd::Zero(rows, cols);
//    T.block(0,0,P.rows(),P.cols()) = P;
//    T.block(P.rows(),P.cols(),A.rows(),A.cols()) = A;
//    T.block(P.rows()+A.rows(),P.cols()+A.cols(),B.rows(),B.cols()) = B;
//}

/** Spin separate this orbital set into subsets of pure spin
 *
 * Separates the orbitals of this set into subsets
 * containing only alpha, beta or paired orbitals.
 * The Fock matrix must be separated accordingly.
 * The ownership of the orbitals is transferred to the
 * new sets and the original set is cleared.
 */
//void OrbitalVector::separateSpinOrbitals(OrbitalVector &phi_p, OrbitalVector &phi_a, OrbitalVector &phi_b) {
//    if (phi_p.size() != this->getNPaired()) MSG_ERROR("Size mismatch paired");
//    if (phi_a.size() != this->getNAlpha()) MSG_ERROR("Size mismatch alpha");
//    if (phi_b.size() != this->getNBeta()) MSG_ERROR("Size mismatch beta");

//    int n_t = this->size();
//    int n_p = 0;
//    int n_a = 0;
//    int n_b = 0;
//    for (int i = 0; i < n_t; i++) {
//        Orbital *orb = getOrbitalPtr(i);
//        if (orb == 0) continue;
//        int spin = orb->getSpin();
//        if (spin == Orbital::Paired) {
//            phi_p.replaceOrbital(n_p, &orb);
//            n_p++;
//        } else if (spin == Orbital::Alpha) {
//            phi_a.replaceOrbital(n_a, &orb);
//            n_a++;
//        } else if (spin == Orbital::Beta) {
//            phi_b.replaceOrbital(n_b, &orb);
//            n_b++;
//        } else {
//            MSG_ERROR("Invalid orbital spin");
//        }
//        this->orbitals[i] = 0;
//    }
//}

//void OrbitalVector::collectSpinOrbitals(OrbitalVector &phi_p, OrbitalVector &phi_a, OrbitalVector &phi_b) {
//    int n_p = phi_p.size();
//    int n_a = phi_a.size();
//    int n_b = phi_b.size();
//    if ((n_p + n_a + n_b) != this->size()) MSG_ERROR("Size mismatch");
//    int n_t = 0;
//    for (int i = 0; i < n_p; i++) {
//        Orbital *orb = phi_p.getOrbitalPtr(i);
//        replaceOrbital(n_t, &orb);
//        n_t++;
//    }
//    for (int i = 0; i < n_a; i++) {
//        Orbital *orb = phi_a.getOrbitalPtr(i);
//        replaceOrbital(n_t, &orb);
//        n_t++;
//    }
//    for (int i = 0; i < n_b; i++) {
//        Orbital *orb = phi_b.getOrbitalPtr(i);
//        replaceOrbital(n_t, &orb);
//        n_t++;
//    }
//    // Do not dealloc orbitals
//    phi_p.clear(false);
//    phi_a.clear(false);
//    phi_b.clear(false);
//}

/** Normalize all orbitals in the set
 */
void OrbitalVector::normalize() {
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb = getOrbital(i);
        orb.normalize();
    }
}

/** Gram-Schmidt orthogonalize the orbitals in the set
 *
 * Nothing happens to the Fock matrix in this process.
 */
//void OrbitalVector::orthogonalize(double prec) {
//    boost::timer rolex;
//    rolex.restart();
//    for (int i = 0; i < this->size(); i++) {
//        Orbital &iOrb = getOrbital(i);
//        for (int j = 0; j < i; j++) {
//            Orbital &jOrb = getOrbital(j);
//            iOrb.orthogonalize(jOrb);
//            iOrb.crop(prec);
//        }
//    }
//    printout(0, "Orthogonalizing                                  ");
//    printout(0, rolex.elapsed() << endl);
//}

/** Orthogonalize all orbitals in this set against all orbitals in the input set
 *
 * Orbitals are NOT orthogonalized within this set
 */
//void OrbitalVector::orthogonalize(double prec, OrbitalVector &orbs) {
//    boost::timer rolex;
//    rolex.restart();
//    for (int i = 0; i < this->size(); i++) {
//        Orbital &iOrb = getOrbital(i);
//        for (int j = 0; j < orbs.size(); j++) {
//            Orbital &jOrb = orbs.getOrbital(j);
//            iOrb.orthogonalize(jOrb);
//            iOrb.crop(prec);
//        }
//    }
//    printout(0, "Orthogonalizing                                  ");
//    printout(0, rolex.elapsed() << endl);
//}

//MatrixXd OrbitalVector::orthonormalize(double prec, MatrixXd *F) {
//    MatrixXd U_p, U_a, U_b;
//    MatrixXd P, A, B;
//    if (F != 0) {
//        separateSpinMatrix(*F, P, A, B);
//    }

//    OrbitalVector phi_p("phi_p", getNDoubly());
//    OrbitalVector phi_a("phi_a", getNAlpha());
//    OrbitalVector phi_b("phi_b", getNBeta());
//    separateSpinOrbitals(phi_p, phi_a, phi_b);

//    if (phi_p.size() > 0) {
//        U_p = phi_p.calcOrthonormalizationMatrix();
//        phi_p.rotate(prec, U_p);
//        if (F != 0) {
//            P = U_p*P*U_p.transpose();
//        }
//    }
//    if (phi_a.size() > 0) {
//        U_a = phi_a.calcOrthonormalizationMatrix();
//        phi_a.rotate(prec, U_a);
//        if (F != 0) {
//            A = U_a*A*U_a.transpose();
//        }
//    }
//    if (phi_b.size() > 0) {
//        U_b = phi_b.calcOrthonormalizationMatrix();
//        phi_b.rotate(prec, U_b);
//        if (F != 0) {
//            B = U_b*B*U_b.transpose();
//        }
//    }
//    if (F != 0) {
//        collectSpinMatrix(*F, P, A, B);
//    }
//    MatrixXd U;
//    collectSpinMatrix(U, U_p, U_a, U_b);
//    collectSpinOrbitals(phi_p, phi_a, phi_b);
//    return U;
//}

//MatrixXd OrbitalVector::localize(double prec, MatrixXd *F) {
//    MatrixXd U_p, U_a, U_b;
//    MatrixXd P, A, B;
//    if (F != 0) {
//        separateSpinMatrix(*F, P, A, B);
//    }

//    OrbitalVector phi_p("phi_p", getNDoubly());
//    OrbitalVector phi_a("phi_a", getNAlpha());
//    OrbitalVector phi_b("phi_b", getNBeta());
//    separateSpinOrbitals(phi_p, phi_a, phi_b);

//    if (phi_p.size() > 0) {
//        U_p = phi_p.calcLocalizationMatrix();
//        phi_p.rotate(prec, U_p);
//        if (F != 0) {
//            P = U_p*P*U_p.transpose();
//        }
//    }
//    if (phi_a.size() > 0) {
//        U_a = phi_a.calcLocalizationMatrix();
//        phi_a.rotate(prec, U_a);
//        if (F != 0) {
//            A = U_a*A*U_a.transpose();
//        }
//    }
//    if (phi_b.size() > 0) {
//        U_b = phi_b.calcLocalizationMatrix();
//        phi_b.rotate(prec, U_b);
//        if (F != 0) {
//            B = U_b*B*U_b.transpose();
//        }
//    }
//    if (F != 0) {
//        collectSpinMatrix(*F, P, A, B);
//    }
//    MatrixXd U;
//    collectSpinMatrix(U, U_p, U_a, U_b);
//    collectSpinOrbitals(phi_p, phi_a, phi_b);
//    return U;
//}

//MatrixXd OrbitalVector::diagonalize(double prec, MatrixXd *F) {
//    if (F == 0) MSG_ERROR("No Fock matrix to diagonalize");
//    boost::timer rolex;
//    rolex.restart();
//    printout(1, "Calculating diagonalization matrix               ");

//    MatrixXd U_p, U_a, U_b;
//    MatrixXd P, A, B;
//    separateSpinMatrix(*F, P, A, B);

//    OrbitalVector phi_p("phi_p", getNPaired());
//    OrbitalVector phi_a("phi_a", getNAlpha());
//    OrbitalVector phi_b("phi_b", getNBeta());
//    separateSpinOrbitals(phi_p, phi_a, phi_b);

//    if (phi_p.size() > 0) {
//        MatrixXd S_tilde = phi_p.calcOverlapMatrix();
//        MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//        P = S_m12*P*S_m12.transpose();

//        MatrixXd M = MathUtils::diagonalizeHermitianMatrix(P);
//        U_p = M.transpose()*S_m12;
//        phi_p.rotate(prec, U_p);
//    }
//    if (phi_a.size() > 0) {
//        MatrixXd S_tilde = phi_a.calcOverlapMatrix();
//        MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//        A = S_m12*A*S_m12.transpose();

//        MatrixXd M = MathUtils::diagonalizeHermitianMatrix(A);
//        U_a = M.transpose()*S_m12;
//        phi_a.rotate(prec, U_a);
//    }
//    if (phi_b.size() > 0) {
//        MatrixXd S_tilde = phi_b.calcOverlapMatrix();
//        MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//        B = S_m12*B*S_m12.transpose();

//        MatrixXd M = MathUtils::diagonalizeHermitianMatrix(B);
//        U_b = M.transpose()*S_m12;
//        phi_b.rotate(prec, U_b);
//    }
//    println(1, rolex.elapsed());

//    MatrixXd U;
//    collectSpinMatrix(*F, P, A, B);
//    collectSpinMatrix(U, U_p, U_a, U_b);
//    collectSpinOrbitals(phi_p, phi_a, phi_b);

//    return U;
//}

/** Calculate overlap matrix within orbital set */
MatrixXcd OrbitalVector::calcOverlapMatrix() {
    OrbitalVector &bra = *this;
    OrbitalVector &ket = *this;
    return bra.calcOverlapMatrix(ket);
}

/** Calculate overlap matrix between two orbital sets */
MatrixXcd OrbitalVector::calcOverlapMatrix(OrbitalVector &ket) {
    OrbitalVector &bra = *this;
    MatrixXcd S = MatrixXcd::Zero(bra.size(), ket.size());
    for (int i = 0; i < bra.size(); i++) {
        Orbital &bra_i = bra.getOrbital(i);
        for (int j = 0; j < ket.size(); j++) {
            Orbital &ket_j = ket.getOrbital(j);
            S(i,j) = bra_i.dot(ket_j);
        }
    }
    return S;
}

/** Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * This operation includes the orthonormalization using the overlap matrix.
 */
//MatrixXd OrbitalVector::calcDiagonalizationMatrix(MatrixXd F) {
//    boost::timer rolex;
//    rolex.restart();
//    printout(1, "Calculating diagonalization matrix               ");

//    MatrixXd S_tilde = calcOverlapMatrix();
//    MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//    F = S_m12*F*S_m12.transpose();

//    MatrixXd M = MathUtils::diagonalizeHermitianMatrix(F);
//    MatrixXd U = M.transpose()*S_m12;

//    println(1, rolex.elapsed());
//    return U;
//}

/** Minimize the spatial extension of orbitals, by a transformation of orbitals
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 *	which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 */
//MatrixXd OrbitalVector::calcLocalizationMatrix() {
//    double optTime = 0.0;
//    double rotTime = 0.0;
//    double totTime = 0.0;

//    mpi::timer totTimer, timer;
//    totTimer.restart();

//    int oldPrec = TelePrompter::setPrecision(5);
//    println(0, endl);

//    timer.restart();
//    RR rr(*this);
//    int nIter = rr.maximize();//compute total U, rotation matrix
//    optTime += timer.elapsed();
//    println(0, "Localized after iteration               " << setw(20) << nIter);

//    MatrixXd U;
//    if (nIter > 0) {
//        U = rr.getTotalU().transpose();
//    } else {
//        timer.restart();
//        U = calcOrthonormalizationMatrix();
//        optTime += timer.elapsed();
//    }

//    println(0, "Calculating rotation matrix                      " << optTime);
//    println(0, "Total time localization                          " << totTimer.elapsed());
//    printout(0, endl);

//    TelePrompter::setPrecision(oldPrec);
//    return U;
//}

//MatrixXd OrbitalVector::calcOrthonormalizationMatrix() {
//    boost::timer rolex;
//    rolex.restart();
//    printout(1, "Calculating orthonormalization matrix            ");

//    MatrixXd S_tilde = calcOverlapMatrix();
//    MatrixXd U = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//    spinCleanMatrix(U);

//    println(1, rolex.elapsed());
//    return U;
//}

//int OrbitalVector::printTreeSizes() const {
//    int nNodes = 0;
//    int nTrees = 0;
//    for (int i = 0; i < size(); i++) {
//        if (this->orbitals[i] != 0) {
//            nNodes += this->orbitals[i]->getNNodes();
//            nTrees++;
//        }
//    }
//    println(0, " OrbitalVector        " << setw(15) << nTrees << setw(25) << nNodes);
//    return nNodes;
//}

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
//RR::RR(OrbitalVector &orbitals) {
//    N = orbitals.size();
//    if (N == 0) {
//        MSG_ERROR("Cannot localize empty set");
//    }
//    total_U = MatrixXd::Identity(N,N);
//    N2h = N*(N-1)/2;
//    gradient = VectorXd(N2h);
//    hessian = MatrixXd(N2h, N2h);
//    r_i_orig = MatrixXd(N,3*N);
//    r_i = MatrixXd(N,3*N);

//    //Make R matrix
//    PositionFunction<3> posFunc;

//    FunctionTree<3> xTree;
//    posFunc.setDir(0);
//    xTree.projectFunction(posFunc);

//    FunctionTree<3> yTree;
//    posFunc.setDir(1);
//    yTree.projectFunction(posFunc);

//    FunctionTree<3> zTree;
//    posFunc.setDir(2);
//    zTree.projectFunction(posFunc);

//    for (int i = 0; i < N; i++) {
//        Orbital &iOrb = orbitals.getOrbital(i);
//        int iSpin = iOrb.getSpin();
//        Orbital xOrb;
//        Orbital yOrb;
//        Orbital zOrb;
//        xOrb.mult(1.0, xTree, 1.0, iOrb, 0);
//        yOrb.mult(1.0, yTree, 1.0, iOrb, 0);
//        zOrb.mult(1.0, zTree, 1.0, iOrb, 0);
//        for (int j = 0; j <= i; j++) {
//            Orbital &jOrb =  orbitals.getOrbital(j);
//            int jSpin = jOrb.getSpin();
//            if (iSpin != jSpin) {
//                MSG_ERROR("Spins must be separated before localization");
//            }
//            r_i_orig(i,j) = jOrb.innerProduct(xOrb);
//            r_i_orig(j,i) = r_i_orig(i,j);
//            r_i_orig(i,j+N) = jOrb.innerProduct(yOrb);
//            r_i_orig(j,i+N) = r_i_orig(i,j+N);
//            r_i_orig(i,j+2*N) = jOrb.innerProduct(zOrb);
//            r_i_orig(j,i+2*N) = r_i_orig(i,j+2*N);
//        }
//    }
//    println(0, "Rotate R matrices");
//    //rotate R matrices into orthonormal basis
//    MatrixXd S_tilde = orbitals.calcOverlapMatrix();
//    MatrixXd S_tilde_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
//    total_U=S_tilde_m12*total_U;
//    MatrixXd r(N, N);
//    for(int dim=0; dim<3; dim++){
//        for (int j=0; j<N; j++) {
//            for (int i=0; i<=j; i++) {
//                r(i,j)=r_i_orig(i,j+dim*N);
//                r(j,i)=r_i_orig(i,j+dim*N);//Enforce symmetry
//            }
//        }
//        r=total_U.transpose()*r*total_U;
//        for (int j=0; j<N; j++) {
//            for (int i=0; i<=j; i++) {
//                r_i(i,j+dim*N)=r(i,j);
//                r_i(j,i+dim*N)=r(i,j);//Enforce symmetry
//            }
//        }
//    }
//    println(0, "done");
//}

/** compute the value of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 */
//double RR::functional() {
//    //s1 is what should be maximized (i.e. the sum of <i R i>^2)
//    double s1=0.0;
//    for (int dim=0; dim<3; dim++) {
//        for (int j=0; j<N; j++) {
//            s1+=r_i(j,j+dim*N)*r_i(j,j+dim*N);
//        }
//    }
//    return s1;
//}

/** Make gradient vector for the RR case
 */
//double RR::make_gradient() {
//    double norm = 0.0;
//    for (int i=0; i<N2h; i++) gradient(i)=0.0 ;
//    for (int dim=0; dim<3; dim++) {
//        int ij=0;
//        for (int j=0; j<N; j++) {
//            for (int i=0; i<j; i++) {
//                gradient(ij)+=4.0*r_i(i,j+dim*N)*(r_i(i,i+dim*N)-r_i(j,j+dim*N));
//                ij++;
//            }
//        }
//    }
//    for (int ij=0; ij<N2h; ij++) {
//        norm += gradient(ij)*gradient(ij);
//    }
//    return sqrt(norm);
//}


/** Make Hessian matrix for the RR case
 */
//double RR::make_hessian() {
//    double djk,djl,dik,dil;
//    for (int j=0; j<N2h; j++){
//        for (int i=0; i<N2h; i++){
//            hessian(i,j)=0.0 ;
//        }
//    }
//    for (int dim=0; dim<3; dim++) {
//        int kl=0;
//        for (int l=0; l<N; l++) {
//            for (int k=0; k<l; k++) {
//                int ij=0;
//                for (int j=0; j<N; j++) {
//                    for (int i=0; i<j; i++) {
//                        djk = j==k ? 1.0: 0.0;
//                        djl = j==l ? 1.0: 0.0;
//                        dik = i==k ? 1.0: 0.0;
//                        dil = i==l ? 1.0: 0.0;

//                        hessian(ij,kl)+=2.0*(
//                                    djk*r_i(i,i+dim*N)*r_i(l,i+dim*N)
//                                    -djl*r_i(i,i+dim*N)*r_i(k,i+dim*N)
//                                    -dik*r_i(j,j+dim*N)*r_i(l,j+dim*N)
//                                    +dil*r_i(j,j+dim*N)*r_i(k,j+dim*N)
//                                    -2*dil*r_i(i,i+dim*N)*r_i(k,j+dim*N)
//                                    +2.0*dik*r_i(i,i+dim*N)*r_i(l,j+dim*N)
//                                    +2*djl*r_i(j,j+dim*N)*r_i(k,i+dim*N)
//                                    -2.0*djk*r_i(j,j+dim*N)*r_i(l,i+dim*N)
//                                    +djk*r_i(l,l+dim*N)*r_i(i,l+dim*N)
//                                    -dik*r_i(l,l+dim*N)*r_i(j,l+dim*N)
//                                    -djl*r_i(k,k+dim*N)*r_i(i,k+dim*N)
//                                    +dil*r_i(k,k+dim*N)*r_i(j,k+dim*N)
//                                    -4*(dil-dik-djl+djk)*r_i(i,j+dim*N)*r_i(k,l+dim*N));

//                        ij++;
//                    }
//                }
//                kl++;
//            }
//        }
//    }

//    return 0;
//}

/** Given the step matrix, update the rotation matrix and the R matrix
 */
//void RR::do_step(VectorXd step){
//    MatrixXd A(N,N);
//    //define rotation U=exp(-A), A real antisymmetric, from step
//    int ij=0;
//    for (int j=0; j<N; j++) {
//        for (int i=0; i<j; i++) {
//            A(i,j)=step(ij);
//            A(j,i)=-A(i,j);
//            ij++;
//        }
//        A(j,j)=0.0;
//    }

//    //calculate U=exp(-A) by diagonalization and U=Vexp(id)Vt with VdVt=iA
//    //could also sum the term in the expansion if A is small
//    total_U*=MathUtils::SkewMatrixExp(A);

//    //rotate the original r matrix with total U
//    MatrixXd r(N, N);
//    for(int dim=0; dim<3; dim++){
//        for (int j=0; j<N; j++) {
//            for (int i=0; i<N; i++) {
//                r(i,j)=r_i_orig(i,j+dim*N);
//            }
//        }
//        r=total_U.transpose()*r*total_U;
//        for (int j=0; j<N; j++) {
//            for (int i=0; i<N; i++) {
//                r_i(i,j+dim*N)=r(i,j);
//            }
//        }
//    }
//}

