//#include <fstream>

#include "Orbital.h"
#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "parallel.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

Orbital* workOrb=0;

Orbital::Orbital(int occ, int s)
        : spin(s),
          occupancy(occ),
          error(1.0),
          real(0),
          imag(0) {
}

Orbital::Orbital(const Orbital &orb)
        : spin(orb.spin),
          occupancy(orb.occupancy),
          error(1.0),
          real(0),
          imag(0) {
}

Orbital& Orbital::operator=(const Orbital &orb) {
    if (this != &orb) {
        if (this->real != 0) MSG_ERROR("Orbital not empty");
        if (this->imag != 0) MSG_ERROR("Orbital not empty");
        this->spin = orb.spin;
        this->occupancy = orb.occupancy;
        this->error = orb.error;
        this->real = orb.real;
        this->imag = orb.imag;
    }
    return *this;
}

void Orbital::clear(bool free) {
    if (this->real != 0 and free) delete this->real;
    if (this->imag != 0 and free) delete this->imag;
    this->real = 0;
    this->imag = 0;
}

void Orbital::allocReal() {
    if (this->real != 0) MSG_ERROR("Orbital not empty");
    this->real = new FunctionTree<3>(*MRA);
}

void Orbital::allocImag() {
    if (this->imag != 0) MSG_ERROR("Orbital not empty");
    this->imag = new FunctionTree<3>(*MRA);
}

int Orbital::getNNodes(int type) const {
    int rNodes = 0;
    int iNodes = 0;
    if (this->hasReal()) rNodes = this->real->getNNodes();
    if (this->hasImag()) iNodes = this->imag->getNNodes();
    if (type == Real) {
        return rNodes;
    }
    if (type == Imag) {
        return iNodes;
    }
    return rNodes + iNodes;
}

double Orbital::getSquareNorm(int type) const {
    double rNorm = 0;
    double iNorm = 0;
    if (this->hasReal()) rNorm = this->real->getSquareNorm();
    if (this->hasImag()) iNorm = this->imag->getSquareNorm();
    if (type == Real) {
        return rNorm;
    }
    if (type == Imag) {
        return iNorm;
    }
    return rNorm + iNorm;
}

void Orbital::compare(const Orbital &orb) const {
    if (this->compareOccupancy(orb) < 0) {
        MSG_WARN("Different occupancy");
    }
    if (this->compareSpin(orb) < 0) {
        MSG_WARN("Different spin");
    }
}

int Orbital::compareOccupancy(const Orbital &orb) const {
    if (this->getOccupancy() == orb.getOccupancy()) {
        return this->getOccupancy();
    }
    return -1;
}

int Orbital::compareSpin(const Orbital &orb) const {
    if (this->getSpin() == orb.getSpin()) {
        return this->getSpin();
    }
    return -1;
}

double Orbital::getExchangeFactor(const Orbital &orb) const {
    if (orb.getSpin() == Paired) {
        return 1.0;
    } else if (this->getSpin() == Paired) {
        return 0.5;
    } else if (this->getSpin() == orb.getSpin()) {
        return 1.0;
    }
    return 0.0;
}

bool Orbital::isConverged(double prec) const {
    if (getError() > prec) {
        return false;
    } else {
        return true;
    }
}

complex<double> Orbital::dot(Orbital &ket) {
    Orbital &bra = *this;
    if ((bra.getSpin() == Alpha) and (ket.getSpin() == Beta)) {
        return 0.0;
    }
    if ((bra.getSpin() == Beta) and (ket.getSpin() == Alpha)) {
        return 0.0;
    }
    double re = 0.0;
    double im = 0.0;
    if (bra.hasReal() and ket.hasReal()) re += bra.re().dot(ket.re());
    if (bra.hasImag() and ket.hasImag()) re += bra.im().dot(ket.im());
    if (bra.hasReal() and ket.hasImag()) im += bra.re().dot(ket.im());
    if (bra.hasImag() and ket.hasReal()) im -= bra.im().dot(ket.re());
    return complex<double>(re, im);
}

void Orbital::normalize() {
    double norm = sqrt(getSquareNorm());
    *this *= 1.0/norm;
}

void Orbital::operator*=(double c) {
    if (hasReal()) (*this->real) *= c;
    if (hasImag()) (*this->imag) *= c;
}

//send an orbital with MPI
void Orbital::send_Orbital(int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  Orbinfo.spin=this->getSpin();
  Orbinfo.occupancy=this->getOccupancy();
  Orbinfo.error=this->getError();
  if(this->hasReal()){
    Orbinfo.NchunksReal = this->re().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksReal = 0;}
  if(this->hasImag()){
    Orbinfo.NchunksImag = this->im().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksImag = 0;}
  
  int count=sizeof(Metadata);
  MPI_Send(&Orbinfo, count, MPI_BYTE, dest, 0, comm);
  
  if(this->hasReal())Send_SerialTree(&this->re(), Orbinfo.NchunksReal, dest, tag, comm);
  if(this->hasImag())Send_SerialTree(&this->im(), Orbinfo.NchunksImag, dest, tag*10000, comm);
  
#endif
}

//send an orbital with MPI
void Orbital::Isend_Orbital(int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  Orbinfo.spin=this->getSpin();
  Orbinfo.occupancy=this->getOccupancy();
  Orbinfo.error=this->getError();
  if(this->hasReal()){
    Orbinfo.NchunksReal = this->re().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksReal = 0;}
  if(this->hasImag()){
    Orbinfo.NchunksImag = this->im().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksImag = 0;}
  
  MPI_Request request;
  int count=sizeof(Metadata);
  MPI_Isend(&Orbinfo, count, MPI_BYTE, dest, 0, comm, &request);
  
  if(this->hasReal())ISend_SerialTree(&this->re(), Orbinfo.NchunksReal, dest, tag, comm);
  if(this->hasImag())ISend_SerialTree(&this->im(), Orbinfo.NchunksImag, dest, tag*10000, comm);
  
#endif
}
//receive an orbital with MPI
void Orbital::Rcv_Orbital(int source, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  int count=sizeof(Metadata);
  MPI_Recv(&Orbinfo, count, MPI_BYTE, source, 0, comm, &status);
  this->setSpin(Orbinfo.spin);
  this->setOccupancy(Orbinfo.occupancy);
  this->setError(Orbinfo.error);
  
  if(Orbinfo.NchunksReal>0){
    if(not this->hasReal()){
      //We must have a tree defined for receiving nodes. Define one:
      this->real = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    Rcv_SerialTree(&this->re(), Orbinfo.NchunksReal, source, tag, comm);}

  if(Orbinfo.NchunksImag>0){
    Rcv_SerialTree(&this->im(), Orbinfo.NchunksImag, source, tag*10000, comm);
  }else{
    //&(this->im())=0;
  }
  
#endif

}

//receive an orbital with MPI
void Orbital::IRcv_Orbital(int source, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;
  MPI_Request request;

  int count=sizeof(Metadata);
  MPI_Irecv(&Orbinfo, count, MPI_BYTE, source, 0, comm, &request);

  this->setSpin(Orbinfo.spin);
  this->setOccupancy(Orbinfo.occupancy);
  this->setError(Orbinfo.error);
  
  if(Orbinfo.NchunksReal>0){
    if(not this->hasReal()){
      //We must have a tree defined for receiving nodes. Define one:
      this->real = new FunctionTree<3>(*MRA,MaxAllocNodes);
    }
    IRcv_SerialTree(&this->re(), Orbinfo.NchunksReal, source, tag, comm);}

  if(Orbinfo.NchunksImag>0){
    IRcv_SerialTree(&this->im(), Orbinfo.NchunksImag, source, tag*10000, comm);
  }else{
    //&(this->im())=0;
  }
  
#endif

}
