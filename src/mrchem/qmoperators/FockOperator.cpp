#include "FockOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "SCFEnergy.h"
#include "MathUtils.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

FockOperator::FockOperator(KineticOperator *t,
                           NuclearPotential *v,
                           CoulombOperator *j,
                           ExchangeOperator *k,
                           XCOperator *xc)
    : QMOperator(MRA->getMaxScale()),
      T(t),
      V(v),
      J(j),
      K(k),
      XC(xc),
      H_1(0) {
}

FockOperator::~FockOperator() {
    this->T = 0;
    this->V = 0;
    this->J = 0;
    this->K = 0;
    this->XC = 0;
    this->H_1 = 0;
}

void FockOperator::setup(double prec) {
    Timer timer;
    TelePrompter::printHeader(0, "Setting up Fock operator");
    if (this->T != 0) this->T->setup(prec);
    if (this->V != 0) this->V->setup(prec);
    if (this->J != 0) this->J->setup(prec);
    if (this->K != 0) this->K->setup(prec);
    if (this->XC != 0) this->XC->setup(prec);
    if (this->H_1 != 0) this->H_1->setup(prec);
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}

void FockOperator::clear() {
    if (this->T != 0) this->T->clear();
    if (this->V != 0) this->V->clear();
    if (this->J != 0) this->J->clear();
    if (this->K != 0) this->K->clear();
    if (this->XC != 0) this->XC->clear();
    if (this->H_1 != 0) this->H_1->clear();
}

void FockOperator::rotate(MatrixXd &U) {
    if (this->K != 0) this->K->rotate(U);
}

Orbital* FockOperator::operator() (Orbital &orb_p) {
    OrbitalAdder add(1.0, this->max_scale);

    vector<Orbital *> orbs;
    if (this->T != 0) orbs.push_back((*this->T)(orb_p));
    if (this->V != 0) orbs.push_back((*this->V)(orb_p));
    if (this->J != 0) orbs.push_back((*this->J)(orb_p));
    if (this->K != 0) orbs.push_back((*this->K)(orb_p));
    if (this->XC != 0) orbs.push_back((*this->XC)(orb_p));
    if (this->H_1 != 0) orbs.push_back((*this->H_1)(orb_p));

    vector<complex<double> > coefs;
    for (int i = 0; i < orbs.size(); i++) coefs.push_back(1.0);

    Orbital *result = new Orbital(orb_p);
    add(*result, coefs, orbs, true);

    for (int n = 0; n < orbs.size(); n++) {
        if (orbs[n] != 0) delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

Orbital* FockOperator::adjoint(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

double FockOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->T != 0) result += (*this->T)(orb_i, orb_j);
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    if (this->H_1 != 0) result += (*this->H_1)(orb_i, orb_j);
    return result;
}

double FockOperator::adjoint(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd FockOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);

#ifdef HAVE_MPI
    for(int i = 0; i<i_orbs.size();i++){
         Orbital& orb_i = i_orbs.getOrbital(i);

	 if(i%MPI_size==MPI_rank){
	   //responsible for this orbital, send it to everybody else. Could use Bcast, but will go another way
	   for(int i_mpi = 0; i_mpi<MPI_size;i_mpi++){
	     if(i_mpi!= MPI_rank)orb_i.send_Orbital(i_mpi, 55);
	   }
	 }else{
	   //get orbital 
	   orb_i.Rcv_Orbital(i%MPI_size, 55);
	 }

	 for(int j = 0; j<j_orbs.size();j++){
	     Orbital &orb_j = j_orbs.getOrbital(j);

	     if(j%MPI_size==MPI_rank){
	       //Only one process does the computations
	       if (this->T != 0) result(i,j) += (*this->T)(orb_i, orb_j);
	       if (this->V != 0) result(i,j) += (*this->V)(orb_i, orb_j);
	       if (this->J != 0) result(i,j) += (*this->J)(orb_i, orb_j);
	       if (this->K != 0) result(i,j) += (*this->K)(orb_i, orb_j);
	       if (this->XC != 0) result(i,j) += (*this->XC)(orb_i, orb_j);
	       if (this->H_1 != 0) result(i,j) += (*this->H_1)(orb_i, orb_j);
	     }
	 }
      }

    MPI_Allreduce(MPI_IN_PLACE, &result(0,0), Ni*Nj,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    
    if (this->T != 0) result += (*this->T)(i_orbs, j_orbs);
    if (this->V != 0) result += (*this->V)(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J)(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K)(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC)(i_orbs, j_orbs);
    if (this->H_1 != 0) result += (*this->H_1)(i_orbs, j_orbs);

#endif
    return result;
}

MatrixXd FockOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);
    if (this->T != 0) result += (*this->T).adjoint(i_orbs, j_orbs);
    if (this->V != 0) result += (*this->V).adjoint(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J).adjoint(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K).adjoint(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC).adjoint(i_orbs, j_orbs);
    if (this->H_1 != 0) result += (*this->H_1).adjoint(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyPotential(Orbital &orb_p) {
    OrbitalAdder add(1.0, this->max_scale);

    vector<Orbital *> orbs;
    if (this->V != 0) orbs.push_back((*this->V)(orb_p));
    if (this->J != 0) orbs.push_back((*this->J)(orb_p));
    if (this->K != 0) orbs.push_back((*this->K)(orb_p));
    if (this->XC != 0) orbs.push_back((*this->XC)(orb_p));

    vector<complex<double> > coefs;
    for (int i = 0; i < orbs.size(); i++) coefs.push_back(1.0);

    Timer timer;
    Orbital *result = new Orbital(orb_p);
    add(*result, coefs, orbs, true);
    timer.stop();
    double time = timer.getWallTime();
    int nNodes = result->getNNodes();
    TelePrompter::printTree(1, "Sum potential operator", nNodes, time);

    for (int n = 0; n < orbs.size(); n++) {
        delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

double FockOperator::applyPotential(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    return result;
}

MatrixXd FockOperator::applyPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    if (this->V != 0) result += (*this->V)(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J)(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K)(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC)(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyAdjointPotential(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

double FockOperator::applyAdjointPotential(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd FockOperator::applyAdjointPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

SCFEnergy FockOperator::trace(OrbitalVector &phi, MatrixXd &F) {
    double E_nuc = 0.0;
    double E_el = 0.0;
    double E_orb = 0.0;
    double E_kin = 0.0;
    double E_en = 0.0;
    double E_ee = 0.0;
    double E_xc = 0.0;
    double E_x = 0.0;

    // Nuclear part
    if (this->V != 0) {
        Nuclei &nucs = this->V->getNuclei();
        int nNucs = nucs.size();
        for (int i = 0; i < nNucs; i++) {
            for (int j = i+1; j < nNucs; j++) {
                const Nucleus &nuc_i = nucs[i];
                const Nucleus &nuc_j = nucs[j];
                double Z_i = nuc_i.getCharge();
                double Z_j = nuc_j.getCharge();
                const double *R_i = nuc_i.getCoord();
                const double *R_j = nuc_j.getCoord();
                double r_ij = MathUtils::calcDistance(3, R_i, R_j);
                E_nuc += (Z_i*Z_j)/r_ij;
            }
        }
    }

    // Electronic part
    double E_xc2 = 0.0;
    if (this->XC != 0) E_xc = this->XC->getEnergy();
    for (int i = 0; i < phi.size(); i++) {
        Orbital &phi_i = phi.getOrbital(i);
        double occ = (double) phi_i.getOccupancy();
        double e_i = occ*F(i,i);
        E_orb += e_i;

        if (this->V != 0) {
            println(2, "\n<" << i << "|V_nuc|" << i << ">");
            E_en += occ*(*this->V)(phi_i,phi_i);
        }
        if (this->J != 0) {
            println(2, "\n<" << i << "|J|" << i << ">");
            E_ee += 0.5*occ*(*this->J)(phi_i,phi_i);
        }
        if (this->K != 0) {
            println(2, "\n<" << i << "|K|" << i << ">");
            E_x += (*this->K)(phi_i,phi_i);
        }
        if (this->XC != 0) {
            println(2, "\n<" << i << "|V_xc|" << i << ">");
            E_xc2 += occ*(*this->XC)(phi_i,phi_i);
        }
    }
    double E_eex = E_ee + E_x;
    double E_orbxc2 = E_orb - E_xc2;
    E_kin = E_orbxc2 - 2.0*E_eex - E_en;
    E_el = E_orbxc2 - E_eex + E_xc;

    return SCFEnergy(E_nuc, E_el, E_orb, E_kin, E_en, E_ee, E_xc, E_x);
}
