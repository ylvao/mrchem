#include "KAIN.h"
#include "Orbital.h"
#include "OrbitalVector.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

/** Calculates the A matrices and b vectors for each orbital individually.
 *
  * \f$ A_{ij} = \langle x^n - x^i | f(x^n) - f(x^j) \rangle \f$
  * \f$ b_{i} = \langle x^n - x^i | f(x^n) \rangle \f$
  *
  * If Fock matrix is included this is treated as an additional ''orbital''
  * and the return vectors have size nOrbs + 1. Frobenius inner product
  * used for the Fock matrix. If separateOrbitals is false the A's and b's
  * are later collected to single entities. */
void KAIN::setupLinearSystem() {
    NOT_IMPLEMENTED_ABORT;
//    int nHistory = this->orbitals.size() - 1;
//    if (nHistory < 1) {
//        MSG_FATAL("Not enough history to setup system of equations");
//    }
//    int nOrbitals = this->orbitals[nHistory]->size();

//    vector<MatrixXd *> A_matrices;
//    vector<VectorXd *> b_vectors;

//    for (int n = 0; n < nOrbitals; n++) {
//        MatrixXd *orbA = new MatrixXd;
//        VectorXd *orbB = new VectorXd;
//        *orbA = MatrixXd::Zero(nHistory, nHistory);
//        *orbB = VectorXd::Zero(nHistory);

//        Orbital &phi_m = this->orbitals[nHistory]->getOrbital(n);
//        Orbital &fPhi_m = this->dOrbitals[nHistory]->getOrbital(n);

//        for (int i = 0; i < nHistory; i++) {
//            Orbital &phi_i = this->orbitals[i]->getOrbital(n);
//            Orbital dPhi_im;
//            dPhi_im.add(1.0, phi_i, -1.0, phi_m);

//            for (int j = 0; j < nHistory; j++) {
//                Orbital &fPhi_j = this->dOrbitals[j]->getOrbital(n);
//                Orbital dfPhi_jm;
//                dfPhi_jm.add(1.0, fPhi_j, -1.0, fPhi_m);

//                // Ref. Harrisons KAIN paper the following has the wrong sign,
//                // but we define the updates (lowercase f) with opposite sign.
//                (*orbA)(i,j) -= dPhi_im.innerProduct(dfPhi_jm);
//            }
//            (*orbB)(i) += dPhi_im.innerProduct(fPhi_m);
//        }
//        A_matrices.push_back(orbA);
//        b_vectors.push_back(orbB);
//    }

//    // Fock matrix is treated as a whole using the Frobenius inner product
//    if (this->orbitals.size() == this->fock.size()) {
//        const MatrixXd &X_m = this->fock[nHistory];
//        const MatrixXd &fX_m = this->dFock[nHistory];

//        MatrixXd *fockA = new MatrixXd;
//        *fockA = MatrixXd::Zero(nHistory, nHistory);

//        VectorXd *fockB = new VectorXd;
//        *fockB = VectorXd::Zero(nHistory);

//        for (int i = 0; i < nHistory; i++) {
//            const MatrixXd &X_i = this->fock[i];
//            MatrixXd dX_im = X_i - X_m;
//            for (int j = 0; j < nHistory; j++) {
//                const MatrixXd &fX_j = this->dFock[j];
//                MatrixXd dfX_jm = fX_j - fX_m;
//                MatrixXd prod = dX_im.transpose()*dfX_jm;
//                (*fockA)(i,j) -= prod.trace();
//            }
//            MatrixXd prod = dX_im.transpose()*fX_m;
//            (*fockB)(i) += prod.trace();
//        }
//        A_matrices.push_back(fockA);
//        b_vectors.push_back(fockB);
//    }

//    sortLinearSystem(A_matrices, b_vectors);

//    for (int n = 0; n < A_matrices.size(); n++) {
//        delete A_matrices[n];
//        delete b_vectors[n];
//    }
}

/** Compute the next step for orbitals and orbital updates given the
  * solution \f$ c \f$ of the linear problem \f$ Ac = b \f$.
  *
  * The next step \f$ \delta x^n \f$ is constructed as
  * \f$ \delta x^n = f(x^n) + \sum_{j=1}^m c_j[(x^j-x^n)+(f(x^j)-f(x^n))]\f$
  */
void KAIN::expandSolution(OrbitalVector &phi,
                          OrbitalVector &dPhi,
                          MatrixXd *F,
                          MatrixXd *dF) {
    NOT_IMPLEMENTED_ABORT;
//    int nHistory = this->orbitals.size() - 1;
//    int nOrbitals = this->orbitals[nHistory]->size();

//    orbs.clear();
//    dOrbs.clear();

//    int m = 0;
//    for (int n = 0; n < nOrbitals; n++) {
//        if (this->sepOrbitals) {
//            m = n;
//        }
//        vector<double> totCoefs;
//        vector<FunctionTree<3> *> totOrbs;

//        Orbital &phi_m = this->orbitals[nHistory]->getOrbital(n);

//        Orbital *orb = new Orbital(phi_m);
//        *orb = phi_m;

//        Orbital &fPhi_m = this->dOrbitals[nHistory]->getOrbital(n);
//        Orbital *orbStep = new Orbital(fPhi_m);
//        totCoefs.push_back(1.0);
//        totOrbs.push_back(&fPhi_m);

//        // Ref. Harrisons KAIN paper the following has the wrong sign,
//        // but we define the updates (lowercase f) with opposite sign
//        // (but not the orbitals themselves).
//        for (int j = 0; j < nHistory; j++) {
//            vector<double> partCoefs;
//            vector<FunctionTree<3> *> partOrbs;
//            Orbital *partStep = new Orbital(*orbStep);

//            Orbital &phi_j = this->orbitals[j]->getOrbital(n);
//            partCoefs.push_back(1.0);
//            partOrbs.push_back(&phi_j);

//            Orbital &fPhi_j = this->dOrbitals[j]->getOrbital(n);
//            partCoefs.push_back(1.0);
//            partOrbs.push_back(&fPhi_j);

//            partCoefs.push_back(-1.0);
//            partOrbs.push_back(&phi_m);
//            partCoefs.push_back(-1.0);
//            partOrbs.push_back(&fPhi_m);
//            partStep->add(partCoefs, partOrbs, 0);

//            double c_j = (*this->c[m])(j);
//            totCoefs.push_back(c_j);
//            totOrbs.push_back(partStep);
//        }

//        orbStep->add(totCoefs, totOrbs, 0);
//        for (int k = 0; k < totOrbs.size(); k++) {
//            // First entry is the last orbital update and should not be deallocated,
//            // all other entries are locally allocated partSteps that must be deleted
//            if (k != 0) delete totOrbs[k];
//            totOrbs[k] = 0;
//        }

//        orbStep->setError(fPhi_m.getError());
//        orbStep->setOccupancy(fPhi_m.getOccupancy());

//        orbs.replaceOrbital(n, &orb);
//        dOrbs.replaceOrbital(n, &orbStep);
//    }

//    // Treat Fock matrix as a whole using Frobenius inner product
//    if (this->fock.size() == this->orbitals.size()) {
//        if (F == 0 or dF == 0) MSG_ERROR("Invalid fock matrix");

//        MatrixXd X_m = this->fock[nHistory];
//        const MatrixXd &fX_m = this->dFock[nHistory];
//        MatrixXd fockStep = MatrixXd::Zero(nOrbitals, nOrbitals);
//        fockStep = fX_m;
//        m = this->c.size();
//        for (int j = 0; j < nHistory; j++) {
//            const MatrixXd &X_j = this->fock[j];
//            const MatrixXd &fX_j = this->dFock[j];
//            MatrixXd tmpX = MatrixXd::Zero(nOrbitals,nOrbitals);
//            tmpX += X_j;
//            tmpX -= X_m;
//            tmpX += fX_j;
//            tmpX -= fX_m;
//            tmpX *= (*this->c[m-1])(j);
//            fockStep += tmpX;
//        }
//        *F = X_m;
//        *dF = fockStep;
//    }
}
