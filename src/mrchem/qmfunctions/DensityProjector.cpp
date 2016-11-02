#include "DensityProjector.h"
#include "OrbitalVector.h"
#include "Density.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;

DensityProjector::DensityProjector(double prec)
    : add(prec, MRA->getMaxScale()),
      mult(prec, MRA->getMaxScale()),
      grid(MRA->getMaxScale()) {
}

void DensityProjector::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void DensityProjector::operator()(Density &rho, Orbital &phi) {
    if (rho.total != 0) MSG_ERROR("Density not empty");
    if (rho.alpha != 0) MSG_ERROR("Density not empty");
    if (rho.beta != 0) MSG_ERROR("Density not empty");

    double occ = 1.0;
    if (not rho.spin) occ = (double) phi.getOccupancy();

    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive multiplication with negative prec");

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        this->grid(*real_2, phi.re());
        this->mult(*real_2, occ, phi.re(), phi.re(), 0);
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        this->grid(*imag_2, phi.im());
        this->mult(*imag_2, occ, phi.im(), phi.im(), 0);
        sum_vec.push_back(imag_2);
    }

    if (rho.spin) {
        if (phi.getSpin() == Paired) {
            rho.alpha = new FunctionTree<3>(*MRA);
            rho.beta = new FunctionTree<3>(*MRA);
            this->grid(*rho.alpha, sum_vec);
            this->grid(*rho.beta, sum_vec);
            this->add(*rho.alpha, sum_vec, 0);
            this->add(*rho.beta, sum_vec, 0);
        }
        if (phi.getSpin() == Alpha) {
            rho.alpha = new FunctionTree<3>(*MRA);
            this->grid(*rho.alpha, sum_vec);
            this->add(*rho.alpha, sum_vec, 0);
            rho.beta = new FunctionTree<3>(*MRA);
            this->grid(*rho.beta, sum_vec);
            rho.beta->setZero();
        }
        if (phi.getSpin() == Beta) {
            rho.beta = new FunctionTree<3>(*MRA);
            this->grid(*rho.beta, sum_vec);
            this->add(*rho.beta, sum_vec, 0);
            rho.alpha = new FunctionTree<3>(*MRA);
            this->grid(*rho.alpha, sum_vec);
            rho.alpha->setZero();
        }
        FunctionTreeVector<3> tot_vec;
        tot_vec.push_back(1.0, rho.alpha);
        tot_vec.push_back(1.0, rho.beta);
        rho.total = new FunctionTree<3>(*MRA);
        this->grid(*rho.total, tot_vec);
        this->add(*rho.total, tot_vec, 0);
    } else {
        rho.total = new FunctionTree<3>(*MRA);
        this->grid(*rho.total, sum_vec);
        this->add(*rho.total, sum_vec, 0);
        rho.alpha = 0;
        rho.beta = 0;
    }
    sum_vec.clear(true);
}

void DensityProjector::operator()(Density &rho, OrbitalVector &phi) {
    if (rho.total != 0) MSG_ERROR("Density not empty");
    if (rho.alpha != 0) MSG_ERROR("Density not empty");
    if (rho.beta != 0) MSG_ERROR("Density not empty");

    int nOrbs = phi.size();
    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive addition with negative prec");
    add.setPrecision(prec/nOrbs);

    FunctionTreeVector<3> total_vec, alpha_vec, beta_vec;
    vector<Density *> dens_vec;
    
    if(MPI_size>1){
      for (int i_Orb = MPI_rank; i_Orb < phi.size(); i_Orb+=MPI_size) {
	Density *rho_i = new Density(rho);	
	if(i_Orb<phi.size()){
	  Orbital &phi_i = phi.getOrbital(i_Orb);
	  (*this)(*rho_i, phi_i);
	  dens_vec.push_back(rho_i);
	  if (rho_i->total != 0) total_vec.push_back(rho_i->total);
	  if (rho_i->alpha != 0) alpha_vec.push_back(rho_i->alpha);
	  if (rho_i->beta != 0) beta_vec.push_back(rho_i->beta);
	}
	for (int iter = 0;  iter<MPI_size ; iter++) {
	  int j_MPI=(MPI_size+iter-MPI_rank)%MPI_size;
	  int j_Orb = j_MPI;
	  Density *rho_j = new Density(rho);
	  if(MPI_rank > j_MPI){
	    //send first bra, then receive ket
	    rho_i->send_Density(j_MPI, i_Orb);
	    rho_j->Rcv_Density(j_MPI, j_Orb);
	  }else if(MPI_rank < j_MPI){
	    //receive first bra, then send ket
	    rho_j->Rcv_Density(j_MPI, j_Orb);
	    rho_i->send_Density(j_MPI, i_Orb);
	  }
	  if(j_MPI%MPI_size!=MPI_rank){
	    dens_vec.push_back(rho_j);
	    if (rho_j->total != 0) total_vec.push_back(rho_j->total);
	    if (rho_j->alpha != 0) alpha_vec.push_back(rho_j->alpha);
	    if (rho_j->beta != 0) beta_vec.push_back(rho_j->beta);
	  }
	}
      }
    }else{
      
      //Serial processing
      for (int i = 0; i < phi.size(); i++) {
	Orbital &phi_i = phi.getOrbital(i);
	Density *rho_i = new Density(rho);
	(*this)(*rho_i, phi_i);
	dens_vec.push_back(rho_i);
	if (rho_i->total != 0) total_vec.push_back(rho_i->total);
	if (rho_i->alpha != 0) alpha_vec.push_back(rho_i->alpha);
	if (rho_i->beta != 0) beta_vec.push_back(rho_i->beta);
      }
    }

    if (not rho.spin) {
        rho.total = new FunctionTree<3>(*MRA);
        if (total_vec.size() > 5) {
            this->add(*rho.total, total_vec);
        } else if (total_vec.size() > 0) {
            this->grid(*rho.total, total_vec);
            this->add(*rho.total, total_vec, 0);
        }
        rho.alpha = 0;
        rho.beta = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }

    for (int i = 0; i < dens_vec.size(); i++) {
        dens_vec[i]->clear();
        delete dens_vec[i];
        dens_vec[i] = 0;
    }

}
