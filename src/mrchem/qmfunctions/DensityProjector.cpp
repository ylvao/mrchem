#include "DensityProjector.h"
#include "OrbitalVector.h"
#include "Density.h"
#include "SerialFunctionTree.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;

DensityProjector::DensityProjector(double prec, int max_scale)
    : add(prec, max_scale),
      mult(prec, max_scale),
      grid(max_scale) {
}

void DensityProjector::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void DensityProjector::operator()(Density &rho, Orbital &phi) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasSpin()) MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta()) MSG_ERROR("Density not empty");

    double occ = 1.0;
    if (not rho.isSpinDensity()) occ = (double) phi.getOccupancy();

    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive multiplication with negative prec");

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        this->grid(*real_2, phi.real());
        this->mult(*real_2, occ, phi.real(), phi.real(), 1);
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        this->grid(*imag_2, phi.imag());
        this->mult(*imag_2, occ, phi.imag(), phi.imag(), 1);
        sum_vec.push_back(imag_2);
    }

    if (rho.isSpinDensity()) {
        if (phi.getSpin() == Orbital::Paired) {
            rho.allocAlpha();
            rho.allocBeta();
            this->grid(rho.alpha(), sum_vec);
            this->grid(rho.beta(), sum_vec);
            this->add(rho.alpha(), sum_vec, 0);
            this->add(rho.beta(), sum_vec, 0);
        }
        if (phi.getSpin() == Orbital::Alpha) {
            rho.allocAlpha();
            this->grid(rho.alpha(), sum_vec);
            this->add(rho.alpha(), sum_vec, 0);

            rho.allocBeta();
            this->grid(rho.beta(), sum_vec);
            rho.beta().setZero();
        }
        if (phi.getSpin() == Orbital::Beta) {
            rho.allocBeta();
            this->grid(rho.beta(), sum_vec);
            this->add(rho.beta(), sum_vec, 0);

            rho.allocAlpha();
            this->grid(rho.alpha(), sum_vec);
            rho.alpha().setZero();
        }
        FunctionTreeVector<3> tot_vec;
        tot_vec.push_back(1.0, &rho.alpha());
        tot_vec.push_back(1.0, &rho.beta());
        rho.allocTotal();
        this->grid(rho.total(), tot_vec);
        this->add(rho.total(), tot_vec, 0);

        FunctionTreeVector<3> spin_vec;
        spin_vec.push_back( 1.0, &rho.alpha());
        spin_vec.push_back(-1.0, &rho.beta());
        rho.allocSpin();
        this->grid(rho.spin(), spin_vec);
        this->add(rho.spin(), spin_vec, 0);
    } else {
        rho.allocTotal();
        this->grid(rho.total(), sum_vec);
        this->add(rho.total(), sum_vec, 0);
    }
    sum_vec.clear(true);
}

void DensityProjector::operator()(Density &rho, OrbitalVector &phi) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta()) MSG_ERROR("Density not empty");

    int nOrbs = phi.size();
    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive addition with negative prec");
    add.setPrecision(prec/nOrbs);

    FunctionTreeVector<3> total_vec, spin_vec, alpha_vec, beta_vec;
    vector<Density *> dens_vec;
    vector<int> rho_i_Ix;
    if(MPI_size>1){
      Density *rho_i;	
      for (int iter = 0;  iter<MPI_size ; iter++) {
	int j_MPI = (MPI_size+iter-MPI_rank)%MPI_size;
	if(MPI_rank > j_MPI){
	  //send first all own bras, then receive all kets from j_MPI
	    int i_Ix = 0;
	  for (int i_Orb = MPI_rank; i_Orb < phi.size(); i_Orb+=MPI_size) {
	    if(iter==0){
	      rho_i = new Density(rho);	
	      Orbital &phi_i = phi.getOrbital(i_Orb);
	      (*this)(*rho_i, phi_i);
	      dens_vec.push_back(rho_i);
	      rho_i_Ix.push_back(dens_vec.size()-1);
	      if (rho_i->hasTotal()) total_vec.push_back(&rho_i->total());
	      if (rho_i->hasSpin()) spin_vec.push_back(&rho_i->spin());
	      if (rho_i->hasAlpha()) alpha_vec.push_back(&rho_i->alpha());
	      if (rho_i->hasBeta()) beta_vec.push_back(&rho_i->beta());
	    }else{
	      rho_i = dens_vec[rho_i_Ix[i_Ix]];
	    }
	    i_Ix++;
	    rho_i->send_Density(j_MPI, i_Orb);
	  }
	  //receive
	  for (int j_Orb = j_MPI; j_Orb < phi.size(); j_Orb+=MPI_size) {
	    Density *rho_j = new Density(rho);
	    rho_j->Rcv_Density(j_MPI, j_Orb);
	    dens_vec.push_back(rho_j);
	    if (rho_j->hasTotal()) total_vec.push_back(&rho_j->total());
	    if (rho_j->hasSpin()) spin_vec.push_back(&rho_j->spin());
	    if (rho_j->hasAlpha()) alpha_vec.push_back(&rho_j->alpha());
	    if (rho_j->hasBeta()) beta_vec.push_back(&rho_j->beta());
	  }
	  
	}else{
	  //receive first all kets from j_MPI then send all own bras
	  if(j_MPI != MPI_rank){
	    for (int j_Orb = j_MPI; j_Orb < phi.size(); j_Orb+=MPI_size) {
	      Density *rho_j = new Density(rho);
	      rho_j->Rcv_Density(j_MPI, j_Orb);
	      if(j_MPI != MPI_rank){
		dens_vec.push_back(rho_j);
	        if (rho_j->hasTotal()) total_vec.push_back(&rho_j->total());
	        if (rho_j->hasSpin()) spin_vec.push_back(&rho_j->spin());
	        if (rho_j->hasAlpha()) alpha_vec.push_back(&rho_j->alpha());
	        if (rho_j->hasBeta()) beta_vec.push_back(&rho_j->beta());
	      }
	    }
	  }
	  //send 
	  int i_Ix = 0;
	  for (int i_Orb = MPI_rank; i_Orb < phi.size(); i_Orb+=MPI_size) {
	    if(iter==0){
	      rho_i = new Density(rho);	
	      Orbital &phi_i = phi.getOrbital(i_Orb);
	      (*this)(*rho_i, phi_i);
	      dens_vec.push_back(rho_i);
	      rho_i_Ix.push_back(dens_vec.size()-1);
	      if (rho_i->hasTotal()) total_vec.push_back(&rho_i->total());
	      if (rho_i->hasSpin()) spin_vec.push_back(&rho_i->spin());
	      if (rho_i->hasAlpha()) alpha_vec.push_back(&rho_i->alpha());
	      if (rho_i->hasBeta()) beta_vec.push_back(&rho_i->beta());
	    }else{
	      rho_i = dens_vec[rho_i_Ix[i_Ix]];
	    }
	    i_Ix++;
	    if(j_MPI != MPI_rank)rho_i->send_Density(j_MPI, i_Orb);
	  }
	}
      }

	/*      for (int i_Orb = MPI_rank; i_Orb < phi.size(); i_Orb+=MPI_size) {
	Density *rho_i = new Density(rho);	
	if(i_Orb<phi.size()){
	  Orbital &phi_i = phi.getOrbital(i_Orb);
	  (*this)(*rho_i, phi_i);
	  dens_vec.push_back(rho_i);
	  if (rho_i->total != 0){ total_vec.push_back(rho_i->total);
		  cout<<MPI_rank<<" push  sze "<<rho_i->total->getSerialFunctionTree()->nodeChunks.size()<<endl;}
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
	}*/
      }else{
      
      //Serial processing
      for (int i = 0; i < phi.size(); i++) {
	Orbital &phi_i = phi.getOrbital(i);
	Density *rho_i = new Density(rho);
	(*this)(*rho_i, phi_i);
	dens_vec.push_back(rho_i);
	if (rho_i->hasTotal()) total_vec.push_back(&rho_i->total());
	if (rho_i->hasSpin()) spin_vec.push_back(&rho_i->spin());
	if (rho_i->hasAlpha()) alpha_vec.push_back(&rho_i->alpha());
	if (rho_i->hasBeta()) beta_vec.push_back(&rho_i->beta());
      }
    } 

    if (total_vec.size() > 5) {
        rho.allocTotal();
        this->add(rho.total(), total_vec);
    } else if (total_vec.size() > 0) {
        rho.allocTotal();
        this->grid(rho.total(), total_vec);
        this->add(rho.total(), total_vec, 0);
    }
    if (spin_vec.size() > 5) {
        rho.allocSpin();
        this->add(rho.spin(), spin_vec);
    } else if (spin_vec.size() > 0) {
        rho.allocSpin();
        this->grid(rho.spin(), spin_vec);
        this->add(rho.spin(), spin_vec, 0);
    }
    if (alpha_vec.size() > 5) {
        rho.allocAlpha();
        this->add(rho.alpha(), alpha_vec);
    } else if (alpha_vec.size() > 0) {
        rho.allocAlpha();
        this->grid(rho.alpha(), alpha_vec);
        this->add(rho.alpha(), alpha_vec, 0);
    }
    if (beta_vec.size() > 5) {
        rho.allocBeta();
        this->add(rho.beta(), beta_vec);
    } else if (beta_vec.size() > 0) {
        rho.allocBeta();
        this->grid(rho.beta(), beta_vec);
        this->add(rho.beta(), beta_vec, 0);
    }

    for (int i = 0; i < dens_vec.size(); i++) {
        dens_vec[i]->clear();
        delete dens_vec[i];
        dens_vec[i] = 0;
    }
}
