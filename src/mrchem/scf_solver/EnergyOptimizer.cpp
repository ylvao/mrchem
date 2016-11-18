#include "EnergyOptimizer.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "HelmholtzOperatorSet.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

extern OrbitalVector workOrbVec;

EnergyOptimizer::EnergyOptimizer(HelmholtzOperatorSet &h)
        : GroundStateSolver(h),
          fOper_np1(0) {
}

EnergyOptimizer::~EnergyOptimizer() {
    if (this->fOper_np1 != 0) MSG_ERROR("Solver not properly cleared");
}

void EnergyOptimizer::setup(FockOperator &fock, OrbitalVector &phi, MatrixXd &F,
                            FockOperator &fock_np1, OrbitalVector &phi_np1) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->fOper_np1 = &fock_np1;

    this->orbitals_n = &phi;
    this->orbitals_np1 = &phi_np1;
    this->dOrbitals_n = new OrbitalVector(phi);
}

void EnergyOptimizer::clear() {
    this->nIter = 0;
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->fOper_np1 = 0;

    delete this->dOrbitals_n;

    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

    resetPrecision();
}

bool EnergyOptimizer::optimize() {
    MatrixXd &F_n = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    bool converged = false;
    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycle();
        adjustPrecision(err_o);

        // Compute electronic energy
        fock.setup(getOrbitalPrecision());
        double E = calcProperty();
        this->property.push_back(E);

        // Iterate Helmholtz operators
        this->helmholtz->initialize(F_n.diagonal());
        applyHelmholtzOperators(phi_np1, F_n, phi_n);
        this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n, true);
        // Compute errors
        VectorXd errors = dPhi_n.getNorms();
#ifdef HAVE_MPI
	//distribute errors among all orbitals
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, &errors(0), 1, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
	
        phi_n.setErrors(errors);
        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);

        // Compute Fock matrix
        MatrixXd F_np1 = F_n + calcFockMatrixUpdate();
        phi_n.clear();
        dPhi_n.clear();
        fock.clear();

        // Rotate orbitals
        if (needLocalization()) {
            localize(fock, F_np1, phi_np1);
        } else if (needDiagonalization()) {
            diagonalize(fock, F_np1, phi_np1);
        } else {
            orthonormalize(fock, F_np1, phi_np1);
        }

        // Update orbitals and Fock matrix
        int nOrbs = phi_np1.size();
        MatrixXd U = MatrixXd::Identity(nOrbs, nOrbs);
        F_n = U.transpose()*F_np1*U;
        this->add.rotate(phi_n, U, phi_np1);
        phi_np1.clear();

        timer.stop();
        printOrbitals(F_n, phi_n);
        printProperty();
        printTimer(timer.getWallTime());

        if (err_p < getPropertyThreshold()) {
            converged = true;
            break;
        }
    }
    printConvergence(converged);
    return converged;
}

MatrixXd EnergyOptimizer::calcFockMatrixUpdate() {
    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");

    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    TelePrompter::printHeader(0,"Computing Fock matrix update");

    Timer timer;
    MatrixXd dS_1;
    MatrixXd dS_2;
    if(MPI_size>1){
      dS_1 = dPhi_n.calcOverlapMatrix_P(phi_n).real();
      dS_2 = phi_np1.calcOverlapMatrix_P(dPhi_n).real();
    }else{
      dS_1 = dPhi_n.calcOverlapMatrix(phi_n).real();
      dS_2 = phi_np1.calcOverlapMatrix(dPhi_n).real();
    }

    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
    XCOperator *xc_n = this->fOper_n->getXCOperator();

    int Ni = phi_np1.size();
    int Nj = dPhi_n.size();
    MatrixXd dV_n = MatrixXd::Zero(Ni,Nj);;
    {   // Nuclear potential matrix is computed explicitly
        Timer timer;
#ifdef HAVE_MPI
	
    int orbVecIx = 0;
    Orbital* orb_i;
    OrbitalVector OrbVecChunk_i(0);
    vector<int> rcv_OrbVec;
    for (int i_Orb = MPI_rank; i_Orb < phi_np1.size(); i_Orb+=MPI_size) {
      for (int iter = 0;  iter<MPI_size ; iter++) {
	int rcv_MPI=(MPI_size+iter-MPI_rank)%MPI_size;
	int rcv_Orb = rcv_MPI+MPI_size*(i_Orb/MPI_size);
	if(MPI_rank > rcv_MPI){
	  //send first bra, then receive ket
	  phi_np1.getOrbital(i_Orb).send_Orbital(rcv_MPI, i_Orb);
	  workOrbVec.getOrbital(orbVecIx).Rcv_Orbital(rcv_MPI, rcv_Orb);
	  orb_i=&workOrbVec.getOrbital(orbVecIx++);
	}else if(MPI_rank < rcv_MPI){
	  //receive first bra, then send ket
	  workOrbVec.getOrbital(orbVecIx).Rcv_Orbital(rcv_MPI, rcv_Orb);
	  phi_np1.getOrbital(i_Orb).send_Orbital(rcv_MPI, i_Orb);
	  orb_i=&workOrbVec.getOrbital(orbVecIx++);
	}else{
	  orb_i=&phi_np1.getOrbital(i_Orb);
	}

	OrbVecChunk_i.push_back(*orb_i);
	rcv_OrbVec.push_back(rcv_Orb);

	if(OrbVecChunk_i.size()>=workOrbVecSize or (iter >= MPI_size-1 and i_Orb+MPI_size>=phi_np1.size())){
	  for (int j = MPI_rank; j < dPhi_n.size(); j+=MPI_size) {
	    
	    Orbital &orb_j = dPhi_n.getOrbital(j);
	    OrbitalVector OrbVecChunk_j(0);
	    OrbVecChunk_j.push_back(orb_j);
	    MatrixXd resultChunk = MatrixXd::Zero(OrbVecChunk_i.size(),OrbVecChunk_j.size());
	    
	    //Only one process does the computations
	    OrbVecChunk_j.getOrbital(0).clear(false);
	    OrbVecChunk_j.getOrbital(0)=orb_j;

	     resultChunk += (*v_n)(OrbVecChunk_i,OrbVecChunk_j);
	    OrbVecChunk_j.clear(false);
	    //copy results into final matrix
	    for (int ix = 0;  ix<OrbVecChunk_i.size() ; ix++) {
	      dV_n(rcv_OrbVec[ix],j) += resultChunk(ix,0);
	    }
	  }
	  OrbVecChunk_i.clearVec(false);
	  rcv_OrbVec.clear();
	  orbVecIx=0;
	}
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &dV_n(0,0), Ni*Nj,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        dV_n = (*v_n)(phi_np1, dPhi_n);
#endif
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Nuclear potential matrix", t);
    }

    MatrixXd F_n;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_n(0, 0, j_n, k_n, xc_n);
        F_n = fock_n(phi_np1, phi_n);
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n", t);
    }

    {   // The n+1 Fock operator needs orthonormalized orbitals
        MatrixXd U = calcOrthonormalizationMatrix(phi_np1);
        this->add.rotate(phi_np1, U);
    }

    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

    println(0,"                                                            ");
    // Do not setup exchange, it must be applied on the fly anyway
    if (j_np1 != 0) j_np1->setup(getOrbitalPrecision());
    if (k_np1 != 0) k_np1->QMOperator::setup(getOrbitalPrecision());
    if (xc_np1 != 0) xc_np1->setup(getOrbitalPrecision());
    println(0,"                                                            ");

    MatrixXd F_np1;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_np1(0, 0, j_np1, k_np1, xc_np1);
        MatrixXd F_1 = fock_np1(phi_n, phi_n);
        MatrixXd F_2 = fock_np1(phi_n, dPhi_n);
        fock_np1.clear();

        F_np1 = F_1 + F_2 + F_2.transpose();
        //MatrixXd F_3 = f_np1(*this->dPhi_n, *this->phi_n);
        //MatrixXd F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
        //MatrixXd F_np1 = F_1 + F_2 + F_3 + F_4;
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n+1", t);
    }

    // Re-computing non-orthogonal phi_np1
    phi_np1.clear();
    this->add(phi_np1, 1.0, phi_n, 1.0, dPhi_n, true);

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd dF_1 = dS_1*(*this->fMat_n);
    MatrixXd dF_2 = dS_2*L;
    MatrixXd dF_3 = F_np1 - F_n;

    // Adding up the pieces
    MatrixXd dF_n = dV_n + dF_1 + dF_2 + dF_3;

    // Symmetrizing
    MatrixXd sym = dF_n + dF_n.transpose();
    dF_n = 0.5 * sym;

    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
    return dF_n;
}
