#include <fstream>

#include "Potential.h"
#include "FunctionTreeVector.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA
extern Orbital* workOrb;

using namespace std;
using namespace Eigen;

Potential::Potential()
    : mult(-1.0),
      real(0),
      imag(0) {
}

Potential::~Potential() {
    if (this->real != 0) MSG_ERROR("Operator not properly deallocated");
    if (this->imag != 0) MSG_ERROR("Operator not properly deallocated");
}

void Potential::setup(double prec) {
    QMOperator::setup(prec);
    this->mult.setPrecision(prec);
}

void Potential::clear() {
    if (this->real != 0) delete this->real;
    if (this->imag != 0) delete this->imag;
    this->real = 0;
    this->imag = 0;
    this->mult.setPrecision(-1.0);
    QMOperator::clear();
}


Orbital* Potential::operator() (Orbital &phi_p) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Timer timer;

    Potential &V = *this;
    Orbital *Vphi_p = new Orbital(phi_p);
    this->mult(*Vphi_p, V, phi_p);

    timer.stop();
    int n = Vphi_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied potential", n, t);
    return Vphi_p;
}

Orbital* Potential::adjoint(Orbital &orb) {
    NOT_IMPLEMENTED_ABORT;
}

double Potential::operator() (Orbital &orb_i, Orbital &orb_j) {
    Orbital *operOrb = (*this)(orb_j);
    complex<double> result = orb_i.dot(*operOrb);
    if (result.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    delete operOrb;
    return result.real();
}

double Potential::adjoint(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd Potential::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    Timer timer;
    TelePrompter::printHeader(1, "Compute Potential Matrix Elements");
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    int ix = 0;
    MatrixXcd M = MatrixXcd::Zero(Ni,Nj);

    if(MPI_size>1){
    VectorXcd MM = VectorXcd::Zero(MPI_size*(Ni*Nj)/MPI_size);//same as M but in another order
    Orbital* orb_j;
    for (int j_Orb = MPI_rank; j_Orb < Nj; j_Orb+=MPI_size) {
        Orbital *operOrb = (*this)(j_orbs.getOrbital(j_Orb));//orbital to send

	//exchange orbitals with others (One at a time)
	for (int iter = 0;  iter<MPI_size ; iter++) {
	  int rcv_MPI=(MPI_size+iter-MPI_rank)%MPI_size;//with who to exchange
	  int rcv_Orb = rcv_MPI+MPI_size*(j_Orb/MPI_size);//witch orbital to receive
	  if(MPI_rank > rcv_MPI){
	    //send first bra, then receive ket
	    operOrb->send_Orbital(rcv_MPI, j_Orb);
	    workOrb->Rcv_Orbital(rcv_MPI, rcv_Orb);
	    orb_j=workOrb;
	  }else if(MPI_rank < rcv_MPI){
	    //receive first bra, then send ket
	    workOrb->Rcv_Orbital(rcv_MPI, rcv_Orb);
	    operOrb->send_Orbital(rcv_MPI, j_Orb);
	    orb_j=workOrb;
	  }else{
	    orb_j=operOrb;
	  }

	  //Project on all own orbitals 
	  for (int i = MPI_rank; i < Ni; i+=MPI_size) {	    
            Orbital &orb_i = i_orbs.getOrbital(i);
            M(i,rcv_Orb) = orb_i.dot(*orb_j);
	    MM((ix++) + MPI_rank*(Ni*Nj)/MPI_size) = M(i,rcv_Orb);
	  }
	}
        delete operOrb;
    }
#ifdef HAVE_MPI
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE_COMPLEX, &MM(0), (Ni*Nj)/MPI_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#endif
    //copy MM into M. do the same loop to get indices right
    for (int I_rank=0; I_rank<MPI_size; I_rank++) {
      ix=0;
      for (int j_Orb = I_rank; j_Orb < Nj; j_Orb+=MPI_size) {
	for (int iter = 0;  iter<MPI_size ; iter++) {
	  int rcv_MPI=(MPI_size+iter-I_rank)%MPI_size;//with who to exchange
	  int rcv_Orb = rcv_MPI+MPI_size*(j_Orb/MPI_size);//wich orbital to receive
	  for (int i = I_rank; i < Ni; i+=MPI_size) {	    
	    M(i,rcv_Orb) = MM((ix++) + I_rank*(Ni*Nj)/MPI_size);
	  }
	}
      }
      }
    //MPI_Allreduce(MPI_IN_PLACE, &M(0,0), Ni*Nj,
    //              MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }else{
      //Serial treatment
      for (int j = 0; j < Nj; j++) {
        Orbital &orb_j = j_orbs.getOrbital(j);
        Orbital *operOrb = (*this)(orb_j);
        for (int i = 0; i < Ni; i++) {
	  Orbital &orb_i = i_orbs.getOrbital(i);
	  M(i,j) = orb_i.dot(*operOrb);
        }
        delete operOrb;
      }
    }

    timer.stop();
    TelePrompter::printFooter(1, timer, 2);
    if (M.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return M.real();
}

MatrixXd Potential::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

int Potential::getNNodes() const {
    int nNodes = 0;
    if (this->real != 0) nNodes += this->real->getNNodes();
    if (this->imag != 0) nNodes += this->imag->getNNodes();
    return nNodes;
}

int Potential::printTreeSizes() const {
    int nNodes = this->getNNodes();
    println(0, " Potential         " << setw(15) << 1 << setw(25) << nNodes);
    return nNodes;
}
