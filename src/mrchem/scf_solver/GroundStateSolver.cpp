#include <Eigen/Eigenvalues>

#include "GroundStateSolver.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "PoissonOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombPotential.h"
#include "ExchangePotential.h"
#include "XCPotential.h"
#include "XCFunctional.h"
#include "SCFEnergy.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "DipoleOperator.h"
#include "MathUtils.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

GroundStateSolver::GroundStateSolver(HelmholtzOperatorSet &h)
        : SCF(h),
          fOper_n(0),
          fMat_n(0),
          orbitals_n(0),
          orbitals_np1(0),
          dOrbitals_n(0) {
}

GroundStateSolver::~GroundStateSolver() {
    if (this->fOper_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->fMat_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->orbitals_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->orbitals_np1 != 0) MSG_ERROR("Solver not properly cleared");
    if (this->dOrbitals_n != 0) MSG_ERROR("Solver not properly cleared");
}

/** Computes the Helmholtz argument for the i-th orbital.
 *
 * Argument contains the potential operator acting on orbital i, and the sum
 * of all orbitals weighted by the Fock matrix. The effect of using inexact
 * Helmholtz operators are included in Lambda, wich is a diagonal matrix
 * with the actual lambda parameters used in the Helmholtz operators.
 *
 * greenArg = \hat{V}orb_i + \sum_j (\Lambda_{ij}-F_{ij})orb_j
 */
Orbital* GroundStateSolver::getHelmholtzArgument(int i,
                                                 MatrixXd &F,
                                                 OrbitalVector &phi,
                                                 bool adjoint) {
    FockOperator &fock = *this->fOper_n;

    double coef = -1.0/(2.0*pi);
    Orbital &phi_i = phi.getOrbital(i);

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd LmF = L - F;

    Orbital *part_1 = fock.applyPotential(phi_i);
    Orbital *part_2;
    if(MPI_size>1){
      part_2 = calcMatrixPart_P(i, LmF, phi);
    }else{
      part_2 = calcMatrixPart(i, LmF, phi);
    }

    if (part_1 == 0) part_1 = new Orbital(phi_i);
    if (part_2 == 0) part_2 = new Orbital(phi_i);

    Timer timer;
    Orbital *arg = new Orbital(phi_i);
    this->add(*arg, coef, *part_1, coef, *part_2, true);

    timer.stop();
    double time = timer.getWallTime();
    int nNodes = arg->getNNodes();
    TelePrompter::printTree(2, "Added arguments", nNodes, time);

    if (part_1 != 0) delete part_1;
    if (part_2 != 0) delete part_2;

    return arg;
}
Orbital* GroundStateSolver::getHelmholtzArgument_1(Orbital &phi_i) {
    Timer timer;
    FockOperator &fock = *this->fOper_n;

    Orbital *part_1 = fock.applyPotential(phi_i);
 
    if (part_1 == 0) part_1 = new Orbital(phi_i);

    timer.stop();
    double time = timer.getWallTime();
    int nNodes = part_1->getNNodes();
    TelePrompter::printTree(2, "Argument 1", nNodes, time);

    return part_1;
}
Orbital* GroundStateSolver::getHelmholtzArgument_2(int i,
						 int*  OrbsIx,
                                                 MatrixXd &F,
                                                 OrbitalVector &phi,
						 Orbital*  part_1,
						 double coef_part1,
						 Orbital &phi_i,
                                                 bool adjoint) {

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd LmF = L - F;

    vector<double> coefs;
    vector<Orbital *> orbs;

    for (int j = 0; j < phi.size(); j++) {
      double coef = LmF(i,OrbsIx[j]);
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

    Orbital *part_2 = new Orbital(phi_i);
    Timer timer;
    if (orbs.size() > 0 ) this->add(*part_2, coefs, orbs, false);

    Orbital *arg = new Orbital(phi_i);
    double coef = -1.0/(2.0*pi);

    this->add(*arg, coef_part1, *part_1, coef, *part_2, true);

    timer.stop();
    double time = timer.getWallTime();
    int nNodes = arg->getNNodes();
    TelePrompter::printTree(2, "Added arguments", nNodes, time);

    if (part_1 != 0) delete part_1;
    if (part_2 != 0) delete part_2;

    return arg;
}

double GroundStateSolver::calcProperty() {
    MatrixXd &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi = *this->orbitals_n;

    SCFEnergy scfEnergy;
    scfEnergy.compute(fock, phi, F);
    this->energy.push_back(scfEnergy);
    return scfEnergy.getElectronicEnergy();
}

double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return fabs(getUpdate(this->property, iter, false));
}

void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();

    TelePrompter::printHeader(0, "                    Energy                 Update      Done ");
    printUpdate(" Kinetic  ",  T_1,  T_1 -  T_0);
    printUpdate(" N-E      ",  V_1,  V_1 -  V_0);
    printUpdate(" Coulomb  ",  J_1,  J_1 -  J_0);
    printUpdate(" Exchange ",  K_1,  K_1 -  K_0);
    printUpdate(" X-C      ", XC_1, XC_1 - XC_0);
    TelePrompter::printSeparator(0, '-');
    printUpdate(" Total    ",  E_1,  E_1 -  E_0);
    TelePrompter::printSeparator(0, '=');
}

/** Minimize the spatial extension of orbitals, by a transformation of orbitals
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 *	which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 */
void GroundStateSolver::localize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    Timer timer;
    RR rr(this->orbPrec[0], phi);
    int n_it = rr.maximize();//compute total U, rotation matrix
    timer.stop();

    MatrixXd U;
    if (nIter > 0) {
        U = rr.getTotalU().transpose();
    } else {
        timer.resume();
        U = calcOrthonormalizationMatrix(phi);
        timer.stop();
    }

    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Calculating localization matrix", n_it, t);
    printout(0, endl);

    F = U*F*U.transpose();
    fock.rotate(U);
    this->add.rotate(phi, U);
}

/** Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * This operation includes the orthonormalization using the overlap matrix.
 */
void GroundStateSolver::diagonalize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    MatrixXd S_m12 = calcOrthonormalizationMatrix(phi);
    F = S_m12.transpose()*F*S_m12;

    Timer timer;
    printout(1, "Calculating diagonalization matrix               ");

    SelfAdjointEigenSolver<MatrixXd> es(F.cols());
    es.compute(F);
    MatrixXd M = es.eigenvectors();
    MatrixXd U = M.transpose()*S_m12;

    timer.stop();
    println(1, timer.getWallTime());

    F = es.eigenvalues().asDiagonal();
    fock.rotate(U);
    this->add.rotate(phi, U);
}

void GroundStateSolver::orthonormalize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    MatrixXd U = calcOrthonormalizationMatrix(phi);
    F = U*F*U.transpose();
    fock.rotate(U);
    this->add.rotate(phi, U);
}

MatrixXd GroundStateSolver::calcOrthonormalizationMatrix(OrbitalVector &phi) {
    Timer timer;
    printout(1, "Calculating orthonormalization matrix            ");

    MatrixXd S_tilde;
    if(MPI_size>1){
      S_tilde = phi.calcOverlapMatrix_P_H(phi).real();
    }else{
      S_tilde = phi.calcOverlapMatrix().real();
    }
    SelfAdjointEigenSolver<MatrixXd> es(S_tilde.cols());
    es.compute(S_tilde);

    MatrixXd A = es.eigenvalues().asDiagonal();
    for (int i = 0; i < A.cols(); i++) {
        A(i,i) = pow(A(i,i), -1.0/2.0);
    }
    MatrixXd B = es.eigenvectors();
    MatrixXd U = B*A*B.transpose();

    timer.stop();
    println(1, timer.getWallTime());
    return U;
}

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
RR::RR(double prec, OrbitalVector &phi) {
    N = phi.size();
    if (N < 2) MSG_ERROR("Cannot localize less than two orbitals");
    total_U = MatrixXd::Identity(N,N);
    N2h = N*(N-1)/2;
    gradient = VectorXd(N2h);
    hessian = MatrixXd(N2h, N2h);
    r_i_orig = MatrixXd(N,3*N);
    r_i = MatrixXd(N,3*N);

    //Make R matrix
    DipoleOperator r_x(0, 0.0);
    DipoleOperator r_y(1, 0.0);
    DipoleOperator r_z(2, 0.0);

    r_x.setup(prec);
    r_y.setup(prec);
    r_z.setup(prec);

    for (int i = 0; i < N; i++) {
        Orbital &phi_i = phi.getOrbital(i);
        int spin_i = phi_i.getSpin();
        for (int j = 0; j <= i; j++) {
            Orbital &phi_j =  phi.getOrbital(j);
            int spin_j = phi_j.getSpin();
            if (spin_i != spin_j) {
                MSG_ERROR("Spins must be separated before localization");
            }
            r_i_orig(i,j) = r_x(phi_i, phi_j);
            r_i_orig(i,j+N) = r_y(phi_i, phi_j);
            r_i_orig(i,j+2*N) = r_z(phi_i, phi_j);
            r_i_orig(j,i) = r_i_orig(i,j);
            r_i_orig(j,i+N) = r_i_orig(i,j+N);
            r_i_orig(j,i+2*N) = r_i_orig(i,j+2*N);
        }
    }
    r_x.clear();
    r_y.clear();
    r_z.clear();

    //rotate R matrices into orthonormal basis
    MatrixXd S_tilde = phi.calcOverlapMatrix().real();
    MatrixXd S_tilde_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);

    total_U=S_tilde_m12*total_U;
    MatrixXd r(N, N);
    for(int dim=0; dim<3; dim++){
        for (int j=0; j<N; j++) {
            for (int i=0; i<=j; i++) {
                r(i,j)=r_i_orig(i,j+dim*N);
                r(j,i)=r_i_orig(i,j+dim*N);//Enforce symmetry
            }
        }
        r=total_U.transpose()*r*total_U;
        for (int j=0; j<N; j++) {
            for (int i=0; i<=j; i++) {
                r_i(i,j+dim*N)=r(i,j);
                r_i(j,i+dim*N)=r(i,j);//Enforce symmetry
            }
        }
    }
}

/** compute the value of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 */
double RR::functional() {
//    //s1 is what should be maximized (i.e. the sum of <i R i>^2)
    double s1=0.0;
    for (int dim=0; dim<3; dim++) {
        for (int j=0; j<N; j++) {
            s1+=r_i(j,j+dim*N)*r_i(j,j+dim*N);
        }
    }
    return s1;
}

/** Make gradient vector for the RR case
 */
double RR::make_gradient() {
    double norm = 0.0;
    for (int i=0; i<N2h; i++) gradient(i)=0.0 ;
    for (int dim=0; dim<3; dim++) {
        int ij=0;
        for (int j=0; j<N; j++) {
            for (int i=0; i<j; i++) {
                gradient(ij)+=4.0*r_i(i,j+dim*N)*(r_i(i,i+dim*N)-r_i(j,j+dim*N));
                ij++;
            }
        }
    }
    for (int ij=0; ij<N2h; ij++) {
        norm += gradient(ij)*gradient(ij);
    }
    return sqrt(norm);
}


/** Make Hessian matrix for the RR case
 */
double RR::make_hessian() {
    double djk,djl,dik,dil;
    for (int j=0; j<N2h; j++){
        for (int i=0; i<N2h; i++){
            hessian(i,j)=0.0 ;
        }
    }
    for (int dim=0; dim<3; dim++) {
        int kl=0;
        for (int l=0; l<N; l++) {
            for (int k=0; k<l; k++) {
                int ij=0;
                for (int j=0; j<N; j++) {
                    for (int i=0; i<j; i++) {
                        djk = j==k ? 1.0: 0.0;
                        djl = j==l ? 1.0: 0.0;
                        dik = i==k ? 1.0: 0.0;
                        dil = i==l ? 1.0: 0.0;

                        hessian(ij,kl)+=2.0*(
                                    djk*r_i(i,i+dim*N)*r_i(l,i+dim*N)
                                    -djl*r_i(i,i+dim*N)*r_i(k,i+dim*N)
                                    -dik*r_i(j,j+dim*N)*r_i(l,j+dim*N)
                                    +dil*r_i(j,j+dim*N)*r_i(k,j+dim*N)
                                    -2*dil*r_i(i,i+dim*N)*r_i(k,j+dim*N)
                                    +2.0*dik*r_i(i,i+dim*N)*r_i(l,j+dim*N)
                                    +2*djl*r_i(j,j+dim*N)*r_i(k,i+dim*N)
                                    -2.0*djk*r_i(j,j+dim*N)*r_i(l,i+dim*N)
                                    +djk*r_i(l,l+dim*N)*r_i(i,l+dim*N)
                                    -dik*r_i(l,l+dim*N)*r_i(j,l+dim*N)
                                    -djl*r_i(k,k+dim*N)*r_i(i,k+dim*N)
                                    +dil*r_i(k,k+dim*N)*r_i(j,k+dim*N)
                                    -4*(dil-dik-djl+djk)*r_i(i,j+dim*N)*r_i(k,l+dim*N));

                        ij++;
                    }
                }
                kl++;
            }
        }
    }

    return 0;
}

/** Given the step matrix, update the rotation matrix and the R matrix
 */
void RR::do_step(VectorXd step){
    MatrixXd A(N,N);
    //define rotation U=exp(-A), A real antisymmetric, from step
    int ij=0;
    for (int j=0; j<N; j++) {
        for (int i=0; i<j; i++) {
            A(i,j)=step(ij);
            A(j,i)=-A(i,j);
            ij++;
        }
        A(j,j)=0.0;
    }

    //calculate U=exp(-A) by diagonalization and U=Vexp(id)Vt with VdVt=iA
    //could also sum the term in the expansion if A is small
    total_U*=MathUtils::SkewMatrixExp(A);

    //rotate the original r matrix with total U
    MatrixXd r(N, N);
    for(int dim=0; dim<3; dim++){
        for (int j=0; j<N; j++) {
            for (int i=0; i<N; i++) {
                r(i,j)=r_i_orig(i,j+dim*N);
            }
        }
        r=total_U.transpose()*r*total_U;
        for (int j=0; j<N; j++) {
            for (int i=0; i<N; i++) {
                r_i(i,j+dim*N)=r(i,j);
            }
        }
    }
}
