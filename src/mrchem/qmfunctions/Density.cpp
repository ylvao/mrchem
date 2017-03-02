#include "Density.h"
#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "TelePrompter.h"
#include "parallel.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

Density::Density(bool s)
    : is_spin(s),
      dens_t(0),
      dens_s(0),
      dens_a(0),
      dens_b(0) {
}

Density::Density(const Density &rho)
    : is_spin(rho.is_spin),
      dens_t(0),
      dens_s(0),
      dens_a(0),
      dens_b(0) {
}

Density::~Density() {
    if (this->hasTotal()) MSG_ERROR("Density not properly deallocated");
    if (this->hasSpin()) MSG_ERROR("Density not properly deallocated");
    if (this->hasAlpha()) MSG_ERROR("Density not properly deallocated");
    if (this->hasBeta()) MSG_ERROR("Density not properly deallocated");
}

void Density::clear() {
    if (this->hasTotal()) delete this->dens_t;
    if (this->hasSpin()) delete this->dens_s;
    if (this->hasAlpha()) delete this->dens_a;
    if (this->hasBeta()) delete this->dens_b;
    this->dens_t = 0;
    this->dens_s = 0;
    this->dens_a = 0;
    this->dens_b = 0;
}

int Density::getNNodes(int type) const {
    int tNodes = 0;
    int sNodes = 0;
    int aNodes = 0;
    int bNodes = 0;
    if (this->hasTotal()) tNodes = this->total().getNNodes();
    if (this->hasSpin()) sNodes = this->spin().getNNodes();
    if (this->hasAlpha()) aNodes = this->alpha().getNNodes();
    if (this->hasBeta()) bNodes = this->beta().getNNodes();
    if (type == Density::Total) return tNodes;
    if (type == Density::Spin) return sNodes;
    if (type == Density::Alpha) return aNodes;
    if (type == Density::Beta) return bNodes;
    return tNodes + sNodes + aNodes + bNodes;
}

void Density::setDensity(int s, FunctionTree<3> *rho) {
    if (s == Density::Total) this->dens_t = rho;
    if (s == Density::Spin) this->dens_s = rho;
    if (s == Density::Alpha) this->dens_a = rho;
    if (s == Density::Beta) this->dens_b = rho;
}

void Density::allocTotal() {
    if (this->hasTotal()) MSG_ERROR("Density not empty");
    this->dens_t = new FunctionTree<3>(*MRA);
}

void Density::allocSpin() {
    if (this->hasSpin()) MSG_ERROR("Density not empty");
    this->dens_s = new FunctionTree<3>(*MRA);
}

void Density::allocAlpha() {
    if (this->hasAlpha()) MSG_ERROR("Density not empty");
    this->dens_a = new FunctionTree<3>(*MRA);
}

void Density::allocBeta() {
    if (this->hasBeta()) MSG_ERROR("Density not empty");
    this->dens_b = new FunctionTree<3>(*MRA);
}

/*
int Density::printTreeSizes() const {
    int nNodes = 0;
    if (this->dens_t != 0) nNodes += this->total().getNNodes();
    if (this->dens_a != 0) nNodes += this->alpha().getNNodes();
    if (this->dens_b != 0) nNodes += this->beta().getNNodes();
    println(0, " Density           " << setw(15) << 1 << setw(25) << nNodes);
    return nNodes;
}
*/


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
  if(this->dens_t){
    Densinfo.NchunksTotal = this->total().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksTotal = 0;}
  if(this->dens_a){
    Densinfo.NchunksAlpha = this->alpha().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksAlpha = 0;}
  if(this->dens_b){
    Densinfo.NchunksBeta = this->beta().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Densinfo.NchunksBeta = 0;}
 

  int count=sizeof(Metadata);
  MPI_Send(&Densinfo, count, MPI_BYTE, dest, tag, comm);
 
  if(this->hasTotal())Send_SerialTree(this->dens_t, Densinfo.NchunksTotal, dest, tag+10000, comm);
  if(this->hasAlpha())Send_SerialTree(this->dens_a, Densinfo.NchunksAlpha, dest, tag+20000, comm);
  if(this->hasBeta())Send_SerialTree(this->dens_b, Densinfo.NchunksBeta, dest, tag+30000, comm);

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
    if(not this->hasTotal()){
      //We must have a tree defined for receiving nodes. Define one:
      this->dens_t = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->dens_t, Densinfo.NchunksTotal, source, tag+10000, comm);}
  if(Densinfo.NchunksAlpha>0){
    if(not this->hasAlpha()){
      //We must have a tree defined for receiving nodes. Define one:
      this->dens_a = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->dens_a, Densinfo.NchunksAlpha, source, tag+20000, comm);}
  if(Densinfo.NchunksBeta>0){
    if(not this->hasBeta()){
      //We must have a tree defined for receiving nodes. Define one:
      this->dens_b = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(this->dens_b, Densinfo.NchunksBeta, source, tag+30000, comm);}
  
#endif
}

