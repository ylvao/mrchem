#include "Density.h"
#include "FunctionTree.h"
#include "TelePrompter.h"
#include "SerialFunctionTree.h"
#include "parallel.h"

using namespace std;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

Density::Density(bool s)
    : spin(s),
      total(0),
      alpha(0),
      beta(0) {
}

Density::Density(const Density &rho)
    : spin(rho.spin),
      total(0),
      alpha(0),
      beta(0) {
}

Density::~Density() {
    if (this->total != 0) MSG_ERROR("Density not properly deallocated");
    if (this->alpha != 0) MSG_ERROR("Density not properly deallocated");
    if (this->beta != 0) MSG_ERROR("Density not properly deallocated");
}

void Density::clear() {
    if (this->total != 0) delete this->total;
    if (this->alpha != 0) delete this->alpha;
    if (this->beta != 0) delete this->beta;
    this->total = 0;
    this->alpha = 0;
    this->beta = 0;
}

int Density::getNNodes() const {
    int nNodes = 0;
    if (this->total != 0) nNodes += this->total->getNNodes();
    if (this->alpha != 0) nNodes += this->alpha->getNNodes();
    if (this->beta != 0) nNodes += this->beta->getNNodes();
    return nNodes;
}

void Density::setDensity(int s, FunctionTree<3> *rho) {
    if (s == Paired) this->total = rho;
    if (s == Alpha) this->alpha = rho;
    if (s == Beta) this->beta = rho;
}

FunctionTree<3>& Density::getDensity(int s) {
    FunctionTree<3> *rho = 0;
    if (s == Paired) rho = this->total;
    if (s == Alpha) rho = this->alpha;
    if (s == Beta) rho = this->beta;
    if (rho == 0) MSG_FATAL("Uninitialized density");
    return *rho;
}

int Density::printTreeSizes() const {
    int nNodes = 0;
    if (this->total != 0) nNodes += this->total->getNNodes();
    if (this->alpha != 0) nNodes += this->alpha->getNNodes();
    if (this->beta != 0) nNodes += this->beta->getNNodes();
    println(0, " Density           " << setw(15) << 1 << setw(25) << nNodes);
    return nNodes;
}


//send Density with MPI
void Density::send_Density(int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    bool spin;
    int NchunksTotal;
    int NchunksAlpha;
    int NchunksBeta;
  };

  Metadata Densinfo;

  Densinfo.spin = this->isSpinDensity();
  if(this->total){
    Densinfo.NchunksTotal = this->total->getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksTotal = 0;}
  if(this->alpha){
    Densinfo.NchunksAlpha = this->alpha->getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksAlpha = 0;}
  if(this->beta){
    Densinfo.NchunksBeta = this->beta->getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksBeta = 0;}
 

  int count=sizeof(Metadata);
  MPI_Send(&Densinfo, count, MPI_BYTE, dest, tag, comm);
 
  if(this->total)Send_SerialTree(this->total, Densinfo.NchunksTotal, dest, tag+10000, comm);
  if(this->alpha)Send_SerialTree(this->alpha, Densinfo.NchunksAlpha, dest, tag+20000, comm);
  if(this->beta)Send_SerialTree(this->beta, Densinfo.NchunksBeta, dest, tag+30000, comm);

#endif
}

//receive Density with MPI
void Density::Rcv_Density(int source, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    bool spin;
    int NchunksTotal;
    int NchunksAlpha;
    int NchunksBeta;
  };

  Metadata Densinfo;

  int count=sizeof(Metadata);
  MPI_Recv(&Densinfo, count, MPI_BYTE, source, tag, comm, &status);

  assert(this->isSpinDensity() == Densinfo.spin);

  if(Densinfo.NchunksTotal>0){
    if(not this->total){
      //We must have a tree defined for receiving nodes. Define one:
      this->total = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->total, Densinfo.NchunksTotal, source, tag+10000, comm);}
  if(Densinfo.NchunksAlpha>0){
    if(not this->alpha){
      //We must have a tree defined for receiving nodes. Define one:
      this->alpha = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->alpha, Densinfo.NchunksAlpha, source, tag+20000, comm);}
  if(Densinfo.NchunksBeta>0){
    if(not this->beta){
      //We must have a tree defined for receiving nodes. Define one:
      this->beta = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->beta, Densinfo.NchunksBeta, source, tag+30000, comm);}
  
#endif

}



//void Density::setup(double prec) {
//    NOT_IMPLEMENTED_ABORT;
//    Timer timer;

//    FunctionTreeVector<3> sq_vec;
//    FunctionTreeVector<3> sum_vec;
//    this->mult.setPrecision(prec);
//    for (int i = 0; i < this->orbitals->size(); i++) {
//        Orbital &phi_i = this->orbitals->getOrbital(i);
//        double occ = (double) phi_i.getOccupancy();
//        if (phi_i.hasReal()) {
//            sq_vec.push_back(phi_i.re());
//            sq_vec.push_back(phi_i.re());
//            FunctionTree<3> *real_2 = this->mult(sq_vec);
//            sq_vec.clear();
//            sum_vec.push_back(occ, *real_2);
//        }
//        if (phi_i.hasImag()) {
//            sq_vec.push_back(phi_i.im());
//            sq_vec.push_back(phi_i.im());
//            FunctionTree<3> *imag_2 = this->mult(sq_vec);
//            sq_vec.clear();
//            sum_vec.push_back(occ, *imag_2);
//        }
//    }
//    if (this->total == 0) {
//        this->add.setPrecision(prec);
//        this->total = this->add(sum_vec);
//    } else {
//        this->add.setPrecision(-1.0);
//        this->clean.setPrecision(prec);
//        int nNodes = this->clean(*this->total);
//        this->add(*this->total, sum_vec);
//    }

//    for (int i = 0; i < sum_vec.size(); i++) {
//        delete sum_vec[i];
//    }
//    sum_vec.clear();

//    timer.stop();
//    int n = this->total->getNNodes();
//    double t = timer.getWallTime();
//    TelePrompter::printTree(1, "Electron density", n, t);
//    TelePrompter::printDouble(1, "Charge integral", this->total->integrate());
//}
