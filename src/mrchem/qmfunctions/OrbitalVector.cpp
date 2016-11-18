#include <Eigen/Eigenvalues>

#include "OrbitalVector.h"
#include "Orbital.h"
#include "parallel.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

extern Orbital workOrb;
OrbitalVector workOrbVec(workOrbVecSize);

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

OrbitalVector& OrbitalVector::operator=(const OrbitalVector &inp) {
    OrbitalVector &out = *this;
    if (&inp != &out) {
        if (inp.size() != out.size()) MSG_ERROR("Size mismatch");
        for (int i = 0; i < out.size(); i++) {
            const Orbital &inp_i = inp.getOrbital(i);
            Orbital &out_i = out.getOrbital(i);
            out_i = inp_i;
        }
    }
    return out;
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
 * Deletes the actual functions in the orbitals and removes them from vector
 * Leaves an empty vector
 */
void OrbitalVector::clearVec(bool free) {
    for (int i = 0; i < this->size(); i++) {
      this->orbitals[i]->clear(free);//must first remove link to trees!
      delete this->orbitals[i];
    }
    this->orbitals.clear();
}

/** Clears each orbital in the vector
 *
 * Deletes the actual functions in the orbitals, keeps
 * the spin and occupancy.
 */
void OrbitalVector::clear(bool free) {
    for (int i = 0; i < this->size(); i++) {
        this->orbitals[i]->clear(free);
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

/** Append an orbital to this set
 *
 */
void OrbitalVector::push_back(Orbital& Orb) {
    Orbital *newOrb = new Orbital(Orb);
    *newOrb = Orb;
    this->orbitals.push_back(newOrb);	
}

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
void OrbitalVector::orthogonalize() {
    NOT_IMPLEMENTED_ABORT;
//    Timer rolex;
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
}

/** Orthogonalize all orbitals in this set against all orbitals in the input set
 *
 * Orbitals are NOT orthogonalized within this set
 */
void OrbitalVector::orthogonalize(OrbitalVector &orbs) {
    NOT_IMPLEMENTED_ABORT;
//    Timer rolex;
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
}

/** Calculate overlap matrix within orbital set */
MatrixXcd OrbitalVector::calcOverlapMatrix() {
    OrbitalVector &bra = *this;
    OrbitalVector &ket = *this;

    //    bra.calcOverlapMatrix_P(ket);//testing
    if(MPI_size>1){
      return bra.calcOverlapMatrix_P_H(ket);
    }else{
      return bra.calcOverlapMatrix(ket);
    }
}

/** Calculate overlap matrix between two orbital sets */
MatrixXcd OrbitalVector::calcOverlapMatrix(OrbitalVector &ket) {
    OrbitalVector &bra = *this;
    Timer tottime;
    MatrixXcd S = MatrixXcd::Zero(bra.size(), ket.size());
    for (int i = 0; i < bra.size(); i++) {
        Orbital &bra_i = bra.getOrbital(i);
        for (int j = 0; j < ket.size(); j++) {
            Orbital &ket_j = ket.getOrbital(j);
            S(i,j) = bra_i.dot(ket_j);
        }
    }
    tottime.stop();
     return S;
}

/** Calculate overlap matrix between two orbital sets using MPI*/
MatrixXcd OrbitalVector::calcOverlapMatrix_P(OrbitalVector &ket) {
#ifdef HAVE_MPI
    OrbitalVector &bra = *this;
    MatrixXcd S = MatrixXcd::Zero(bra.size(), ket.size());
    MatrixXcd S_MPI = MatrixXcd::Zero(bra.size(), ket.size());
    assert(bra.size()==ket.size());
    int N = bra.size();
    int NiterMax = ((N+MPI_size-1)/MPI_size)*((N+MPI_size-1)/MPI_size)*MPI_size;
    int doi[NiterMax];
    int doj[NiterMax];
    int sendto[NiterMax];
    int sendorb[NiterMax];
    int rcvfrom[NiterMax];
    int rcvorb[NiterMax];
    int MaxIter=NiterMax;
    Assign_NxN(N, doi, doj, sendto, sendorb, rcvorb, &MaxIter);
 
    Orbital* bra_i;//NB: empty
    Orbital* ket_j;//NB: empty
    Orbital* ket_i;//NB: empty
    Orbital* bra_j;//NB: empty
    //    Orbital rcvOrb(bra.getOrbital(0));//NB: empty
    Timer timer, timerw, tottime;
    timerw.stop();
    timer.stop();

    for (int iter = 0; iter <= MaxIter; iter++) {
    timer.resume();
      int i = doi[iter];
      int j = doj[iter];
      if(sendorb[iter]>=0){
	assert(sendto[iter]>=0);
	if(i%MPI_size >= j%MPI_size){
	  bra_i = &(bra.getOrbital(i));
	  //send first bra, then receive 
	  bra_i->send_Orbital(sendto[iter], sendorb[iter]);
	  if(rcvorb[iter]>=0)workOrb.Rcv_Orbital(rcvorb[iter]%MPI_size, rcvorb[iter]);
	}else{
	  bra_i = &(bra.getOrbital(i));
	  //receive first bra, then send 
	  if(rcvorb[iter]>=0)workOrb.Rcv_Orbital(rcvorb[iter]%MPI_size, rcvorb[iter]);
	  bra_i->send_Orbital(sendto[iter], sendorb[iter]);
	}
      }else if(rcvorb[iter]>=0){
	//receive only, do not send anything
	workOrb.Rcv_Orbital(rcvorb[iter]%MPI_size, rcvorb[iter]);
      }
    timer.stop();
    timerw.resume();
	  
      //sendings finished, do the work
      if(doi[iter]>=0 ){
	assert(i%MPI_size==MPI_rank);//in present implementation, i is always owned locally
	assert(doj[iter]>=0);
	ket_i = &(ket.getOrbital(i));
	if(j%MPI_size!=MPI_rank){bra_j=&workOrb;
	  }else{bra_j= &(bra.getOrbital(j));}
            S_MPI(j,i) =  bra_j->dot(*ket_i);	
      }
      timerw.stop();
   }

    MPI_Allreduce(MPI_IN_PLACE, &S_MPI(0,0), N*N,
                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return S_MPI;
#else
    NOT_REACHED_ABORT;
#endif
}
/** Calculate overlap matrix between two orbital sets using MPI
 * assumes Hermitian overlap
 */
MatrixXcd OrbitalVector::calcOverlapMatrix_P_H(OrbitalVector &ket) {
#ifdef HAVE_MPI
    OrbitalVector &bra = *this;
    MatrixXcd S = MatrixXcd::Zero(bra.size(), ket.size());
    MatrixXcd S_MPI = MatrixXcd::Zero(bra.size(), ket.size());
    assert(bra.size()==ket.size());
    int N = bra.size();
    int NiterMax = ((N+MPI_size-1)/MPI_size)*((N+MPI_size-1)/MPI_size)*MPI_size;
    int doi[NiterMax];
    int doj[NiterMax];
    int sendto[NiterMax];
    int sendorb[NiterMax];
    int rcvfrom[NiterMax];
    int rcvorb[NiterMax];
    int MaxIter=NiterMax;
    Assign_NxN_sym(N, doi, doj, sendto, sendorb, rcvorb, &MaxIter);
 
    Orbital* myOrb_i;//NB: empty
    Orbital* Orb_j;//NB: empty

    //  Orbital rcvOrb(bra.getOrbital(0));//NB: empty
    Timer timer, timerw, tottime;
    timerw.stop();
    timer.stop();
    
    for (int iter = 0; iter <= MaxIter; iter++) {
      timer.resume();
      int i = doi[iter];
      int j = doj[iter];
      //send first , then receive 
      if(sendorb[iter]>=0){
	assert(sendto[iter]>=0);
	myOrb_i = &(bra.getOrbital(i));
	myOrb_i->Isend_Orbital(sendto[iter], sendorb[iter]);
      }
      if(rcvorb[iter]>=0){
	//receive 
	workOrb.Rcv_Orbital(rcvorb[iter]%MPI_size, rcvorb[iter]);
      }
      timer.stop();
      timerw.resume();
      
      //sendings finished, do the work
      if(doi[iter]>=0 ){
	assert(i%MPI_size==MPI_rank);//in present implementation, i is always owned locally
	assert(doj[iter]>=0);
	myOrb_i = &(bra.getOrbital(i));
	if(j%MPI_size!=MPI_rank){
	  Orb_j = &workOrb;
	}else{
	  Orb_j= &(ket.getOrbital(j));
	}
	S_MPI(i,j) =  myOrb_i->dot(*Orb_j);	
	S_MPI(j,i) =  conj(S_MPI(i,j));	
     }
      timerw.stop();
    }
    Timer t1;
    MPI_Allreduce(MPI_IN_PLACE, &S_MPI(0,0), N*N,
                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    t1.stop();
    tottime.stop();
    //if(MPI_rank==0)cout<<" time orbital send/rcv "<<timer<<" MPI_reduce "<<t1<<" overlap "<<timerw<<" total "<< tottime<<endl;

    return S_MPI;
#else
    NOT_REACHED_ABORT;
#endif
}

int OrbitalVector::printTreeSizes() const {
    int nNodes = 0;
    int nTrees = 0;
    for (int i = 0; i < size(); i++) {
        if (this->orbitals[i] != 0) {
            nNodes += this->orbitals[i]->getNNodes();
            nTrees++;
        }
    }
    println(0, " OrbitalVector     " << setw(15) << nTrees << setw(25) << nNodes);
    return nNodes;
}

