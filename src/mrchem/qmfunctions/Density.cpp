#include "Density.h"
#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "ProjectedNode.h"
#include "TelePrompter.h"
#include "parallel.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

Density::Density(bool spin, bool shared)
    : is_spin(spin),
      is_shared(shared),
      dens_t(0),
      dens_s(0),
      dens_a(0),
      dens_b(0) {
    if(this->is_shared and MPI_SH_size>1){
	//Ideally shared memory should be allocated here
    }else{
	this->setIsShared(false);
    }
}
Density::Density(bool s)
    : is_spin(s),
      is_shared(false),
      dens_t(0),
      dens_s(0),
      dens_a(0),
      dens_b(0) {
}

Density::Density(const Density &rho)
    : is_spin(rho.is_spin),
      is_shared(false),//shared memory is not transfered, but must be set explicitely
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
    if(this->isShared()){
	//must first free the shared memory
	MPI_Win_free(&this->MPI_win);//NB: collective command
    }
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
    this->dens_t->getSerialFunctionTree()->isShared = this->isShared();
}

void Density::allocSpin() {
    if (this->hasSpin()) MSG_ERROR("Density not empty");
    this->dens_s = new FunctionTree<3>(*MRA);
    this->dens_s->getSerialFunctionTree()->isShared = this->isShared();
}

void Density::allocAlpha() {
    if (this->hasAlpha()) MSG_ERROR("Density not empty");
    this->dens_a = new FunctionTree<3>(*MRA);
    this->dens_a->getSerialFunctionTree()->isShared = this->isShared();
}

void Density::allocBeta() {
    if (this->hasBeta()) MSG_ERROR("Density not empty");
    this->dens_b = new FunctionTree<3>(*MRA);
    this->dens_b->getSerialFunctionTree()->isShared = this->isShared();
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
    MPI_Comm comm=MPI_Comm_Orb;

    struct Metadata{
	bool spin;
	int NchunksTotal;
	int NchunksSpin;
	int NchunksAlpha;
	int NchunksBeta;
    };

    Metadata Densinfo;

    Densinfo.spin = this->isSpinDensity();
    if(this->dens_t){
	Densinfo.NchunksTotal = this->total().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
    }else{Densinfo.NchunksTotal = 0;}
    if(this->dens_s){
	Densinfo.NchunksSpin = this->spin().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
    }else{Densinfo.NchunksSpin = 0;}
    if(this->dens_a){
	Densinfo.NchunksAlpha = this->alpha().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
    }else{Densinfo.NchunksAlpha = 0;}
    if(this->dens_b){
	Densinfo.NchunksBeta = this->beta().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
    }else{Densinfo.NchunksBeta = 0;}
 

    int count=sizeof(Metadata);
    MPI_Send(&Densinfo, count, MPI_BYTE, dest, tag, comm);
    if(this->hasTotal())Send_SerialTree(this->dens_t, Densinfo.NchunksTotal, dest, tag+10000, comm);
    if(this->hasSpin())Send_SerialTree(this->dens_s, Densinfo.NchunksSpin, dest, tag+15000, comm);
    if(this->hasAlpha())Send_SerialTree(this->dens_a, Densinfo.NchunksAlpha, dest, tag+20000, comm);
    if(this->hasBeta())Send_SerialTree(this->dens_b, Densinfo.NchunksBeta, dest, tag+30000, comm);

#endif
}

//receive Density with MPI
void Density::Rcv_Density(int source, int tag){
#ifdef HAVE_MPI
    MPI_Status status;
    MPI_Comm comm=MPI_Comm_Orb;

    struct Metadata{
	bool spin;
	int NchunksTotal;
	int NchunksSpin;
	int NchunksAlpha;
	int NchunksBeta;
    };

    Metadata Densinfo;

    int count=sizeof(Metadata);
    MPI_Recv(&Densinfo, count, MPI_BYTE, source, tag, comm, &status);

    if(Densinfo.NchunksTotal>0){
	if (not this->hasTotal()) this->allocTotal();
	Rcv_SerialTree(this->dens_t, Densinfo.NchunksTotal, source, tag+10000, comm);}
    if(Densinfo.NchunksSpin>0){
	if (not this->hasSpin()) this->allocSpin();
	Rcv_SerialTree(this->dens_s, Densinfo.NchunksSpin, source, tag+15000, comm);}
    if(Densinfo.NchunksAlpha>0){
	if (not this->hasAlpha()) this->allocAlpha();
	Rcv_SerialTree(this->dens_a, Densinfo.NchunksAlpha, source, tag+20000, comm);}
    if(Densinfo.NchunksBeta>0){
	if (not this->hasBeta()) this->allocBeta();
	Rcv_SerialTree(this->dens_b, Densinfo.NchunksBeta, source, tag+30000, comm);}
  
#endif
}
//Allocate an empty Density where the coefficients are shared between processes within a compute-node
void Density::Allocate_Shared_Density(int shared_size){
    //NB: for now only total spin is shared
#ifdef HAVE_MPI
    double * shared_mem_ptr;
    //allocate a big chunk of shared memory
    Share_memory(shared_size, shared_mem_ptr, this->MPI_win);

    //  if(this->hasTotal()){
    SerialFunctionTree<3>* serialTree_p = (SerialFunctionTree<3>*) this->dens_t->getSerialTree();
    if(serialTree_p->nodeCoeffChunks.size()!=1)MSG_ERROR("shared density must start from empty density");
    
    //repoint the root coeff
    MWNode<3> **roots = this->dens_t->getRootBox().getNodes();
    int nRoots = this->dens_t->getRootBox().size();
    double * d_ptr = shared_mem_ptr;
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
	roots[rIdx]->coefs =  d_ptr; //repoint adress of root nodecoeff
	d_ptr += serialTree_p->sizeNodeCoeff; 
    }
    //deallocate the old coefficients
    delete[] serialTree_p->nodeCoeffChunks[0];//deallocate old chunk
    
    //preallocate all chunks
    int size_NodeCoeffChunk = serialTree_p->sizeNodeCoeff*serialTree_p->maxNodesPerChunk;//in units of doubles
    int size_NodeChunk = sizeof(ProjectedNode<3>)*serialTree_p->maxNodesPerChunk;//in Bytes
    d_ptr = shared_mem_ptr;
    serialTree_p->nodeCoeffChunks[0]= shared_mem_ptr;//reset adress of first chunk
   
    long int nextsize = 2*size_NodeCoeffChunk;//already set one, and need space for next.
    long int availsize = shared_size;
    availsize *= 1024*1024;
    int maxNChunks = serialTree_p->maxNodes / serialTree_p->maxNodesPerChunk -1;

    for (int i = 0; nextsize < availsize and i  <maxNChunks; i++) {
	d_ptr += size_NodeCoeffChunk;
	nextsize += size_NodeCoeffChunk;
	serialTree_p->nodeCoeffChunks.push_back(d_ptr);
	serialTree_p->sNodes = (ProjectedNode<3>*) new char[size_NodeChunk];//NB: nodes (without coeff) are allocated individually, not shared
	serialTree_p->nodeChunks.push_back(serialTree_p->sNodes);
    }
    // }
    /*  if(this->hasSpin()){
	SerialFunctionTree<3>* serialTree_p = (SerialFunctionTree<3>*) this->dens_s->getSerialTree();
	if(serialTree_p->nodeCoeffChunks.size()!=1)MSG_ERROR("shared density must start from empty density");
    
	//repoint the root coeff
	MWNode<3> **roots = this->dens_s->getRootBox().getNodes();
	int nRoots = this->dens_s->getRootBox().size();
	double * d_ptr = shared_mem_ptr;
	for (int rIdx = 0; rIdx < nRoots; rIdx++) {
	roots[rIdx]->coefs =  d_ptr; //repoint adress of root nodecoeff
	d_ptr += serialTree_p->sizeNodeCoeff; 
	}
	//deallocate the old coefficients
	delete[] serialTree_p->nodeCoeffChunks[0];//deallocate old chunk
    
	//preallocate all chunks
	int size_NodeCoeffChunk = serialTree_p->sizeNodeCoeff*serialTree_p->maxNodesPerChunk;//in units of doubles
	int size_NodeChunk = sizeof(ProjectedNode<3>)*serialTree_p->maxNodesPerChunk;//in Bytes
	d_ptr = shared_mem_ptr;
	serialTree_p->nodeCoeffChunks[0]= shared_mem_ptr;//reset adress of first chunk
    
	long int nextsize = 2*size_NodeCoeffChunk;//already set one, and need space for next.
	long int availsize = shared_size;
	availsize *= 1024*1024;
	int maxNChunks = serialTree_p->maxNodes / serialTree_p->maxNodesPerChunk -1;
	for (int i = 0; nextsize < availsize and i  <maxNChunks; i++) {
	d_ptr += size_NodeCoeffChunk;
	nextsize += size_NodeCoeffChunk;
	serialTree_p->nodeCoeffChunks.push_back(d_ptr);
	serialTree_p->sNodes = (ProjectedNode<3>*) new char[size_NodeChunk];//NB: nodes (without coeff) are allocated individually, not shared
	serialTree_p->nodeChunks.push_back(serialTree_p->sNodes);
	}
	}

	if(this->hasAlpha()){
	SerialFunctionTree<3>* serialTree_p = (SerialFunctionTree<3>*) this->dens_a->getSerialTree();
	if(serialTree_p->nodeCoeffChunks.size()!=1)MSG_ERROR("shared density must start from empty density");
    
	//repoint the root coeff
	MWNode<3> **roots = this->dens_a->getRootBox().getNodes();
	int nRoots = this->dens_a->getRootBox().size();
	double * d_ptr = shared_mem_ptr;
	for (int rIdx = 0; rIdx < nRoots; rIdx++) {
	roots[rIdx]->coefs =  d_ptr; //repoint adress of root nodecoeff
	d_ptr += serialTree_p->sizeNodeCoeff; 
	}
	//deallocate the old coefficients
	delete[] serialTree_p->nodeCoeffChunks[0];//deallocate old chunk
    
	//preallocate all chunks
	int size_NodeCoeffChunk = serialTree_p->sizeNodeCoeff*serialTree_p->maxNodesPerChunk;//in units of doubles
	int size_NodeChunk = sizeof(ProjectedNode<3>)*serialTree_p->maxNodesPerChunk;//in Bytes
	d_ptr = shared_mem_ptr;
	serialTree_p->nodeCoeffChunks[0]= shared_mem_ptr;//reset adress of first chunk
    
	long int nextsize = 2*size_NodeCoeffChunk;//already set one, and need space for next.
	long int availsize = shared_size;
	availsize *= 1024*1024;
	int maxNChunks = serialTree_p->maxNodes / serialTree_p->maxNodesPerChunk -1;
	for (int i = 0; nextsize < availsize and i  <maxNChunks; i++) {
	d_ptr += size_NodeCoeffChunk;
	nextsize += size_NodeCoeffChunk;
	serialTree_p->nodeCoeffChunks.push_back(d_ptr);
	serialTree_p->sNodes = (ProjectedNode<3>*) new char[size_NodeChunk];//NB: nodes (without coeff) are allocated individually, not shared
	serialTree_p->nodeChunks.push_back(serialTree_p->sNodes);
	}
	}
	if(this->hasBeta()){
	SerialFunctionTree<3>* serialTree_p = (SerialFunctionTree<3>*) this->dens_b->getSerialTree();
	if(serialTree_p->nodeCoeffChunks.size()!=1)MSG_ERROR("shared density must start from empty density");
    
	//repoint the root coeff
	MWNode<3> **roots = this->dens_b->getRootBox().getNodes();
	int nRoots = this->dens_b->getRootBox().size();
	double * d_ptr = shared_mem_ptr;
	for (int rIdx = 0; rIdx < nRoots; rIdx++) {
	roots[rIdx]->coefs =  d_ptr; //repoint adress of root nodecoeff
	d_ptr += serialTree_p->sizeNodeCoeff; 
	}
	//deallocate the old coefficients
	delete[] serialTree_p->nodeCoeffChunks[0];//deallocate old chunk
    
	//preallocate all chunks
	int size_NodeCoeffChunk = serialTree_p->sizeNodeCoeff*serialTree_p->maxNodesPerChunk;//in units of doubles
	int size_NodeChunk = sizeof(ProjectedNode<3>)*serialTree_p->maxNodesPerChunk;//in Bytes
	d_ptr = shared_mem_ptr;
	serialTree_p->nodeCoeffChunks[0]= shared_mem_ptr;//reset adress of first chunk
    
	long int nextsize = 2*size_NodeCoeffChunk;//already set one, and need space for next.
	long int availsize = shared_size;
	availsize *= 1024*1024;
	int maxNChunks = serialTree_p->maxNodes / serialTree_p->maxNodesPerChunk -1;
	for (int i = 0; nextsize < availsize and i  <maxNChunks; i++) {
	d_ptr += size_NodeCoeffChunk;
	nextsize += size_NodeCoeffChunk;
	serialTree_p->nodeCoeffChunks.push_back(d_ptr);
	serialTree_p->sNodes = (ProjectedNode<3>*) new char[size_NodeChunk];//NB: nodes (without coeff) are allocated individually, not shared
	serialTree_p->nodeChunks.push_back(serialTree_p->sNodes);
	}
	}*/
#else
    MSG_ERROR("Cannot allocate shared memory for serial runs");
#endif
}


