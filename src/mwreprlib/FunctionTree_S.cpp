/**
 *  \date Oct 12, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 */

#include "MWTree_S.h"
#include "FunctionTree_S.h"
#include "MRNode.h"
#include "MWNode.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
//#include "HilbertIterator.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
template<int D>
FunctionTree_S<D>::FunctionTree_S(const MultiResolutionAnalysis<D> &mra, int MaxNumberofNodes)
  : MWTree_S<D> (mra, MaxNumberofNodes) {
    this->LastNode = this->AllocNodes(this->MWTree_p->getRootBox().size());
    for (int rIdx = 0; rIdx < this->MWTree_p->getRootBox().size(); rIdx++) {
        const NodeIndex<D> &nIdx = this->MWTree_p->getRootBox().getNodeIndex(rIdx);
	//        MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
        MRNode<D> *  f_Node  = new (this->LastNode) ProjectedNode<D>(*((FunctionTree<D>*) this->MWTree_p), nIdx);
        this->MWTree_p->getRootBox().setNode(rIdx, &f_Node);
	this->LastNode++;
    }
    this->MWTree_p->resetEndNodeTable();
    this->MWTree_p->calcSquareNorm();
}

/** FunctionTree copy constructor.
  * Copy polynomial order and type, as well as the world box from the
  * given tree, but only at root scale. Initializes the function to zero.
  * Use = operator to copy data.*/
template<int D>
FunctionTree_S<D>::FunctionTree_S(const MWTree_S<D> &tree)
        : MWTree_S<D> (tree) {
    this->LastNode = this->AllocNodes(this->MWTree_p->getRootBox().size());
    for (int rIdx = 0; rIdx < this->MWTree_p->getRootBox().size(); rIdx++) {
        const NodeIndex<D> &nIdx = this->MWTree_p->getRootBox().getNodeIndex(rIdx);
        //MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
	MRNode<D> *  f_tree_S_start  = new (this->LastNode) ProjectedNode<D>(*((FunctionTree<D>*) this->MWTree_p), nIdx);
       this->MWTree_p->getRootBox().setNode(rIdx, &f_tree_S_start);
       this->LastNode++;
    }
    this->MWTree_p->resetEndNodeTable();
    this->MWTree_p->calcSquareNorm();
}

/** FunctionTree copy constructor.
  * Copy polynomial order and type, as well as the world box from the
  * given tree, but only at root scale. Initializes the function to zero.
  * Use = operator to copy data.*/
template<int D>
FunctionTree_S<D>::FunctionTree_S(const FunctionTree_S<D> &tree)
        : MWTree_S<D> (tree) {
    this->LastNode = this->AllocNodes(this->MWTree_p->getRootBox().size());
    for (int rIdx = 0; rIdx < this->MWTree_p->getRootBox().size(); rIdx++) {
      const NodeIndex<D> &nIdx = this->MWTree_p->getRootBox().getNodeIndex(rIdx);
	//        MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
        MRNode<D> *  f_tree_S_start  = new (this->LastNode) ProjectedNode<D>(*((FunctionTree<D>*) this->MWTree_p), nIdx);
        this->MWTree_p->getRootBox().setNode(rIdx, &f_tree_S_start);
	this->LastNode++;
    }
    this->MWTree_p->resetEndNodeTable();
    this->MWTree_p->calcSquareNorm();
}

template<int D>
FunctionTree_S<D>& FunctionTree_S<D>::operator=(const FunctionTree_S<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

/** FunctionTree destructor. Not relevant for MWTree_S*/
template<int D>
FunctionTree_S<D>::~FunctionTree_S() {
  
}

/** Leaves the tree inn the same state as after construction*/
template<int D>
void FunctionTree_S<D>::clear() {
    NOT_IMPLEMENTED_ABORT;
}

/** Write the tree structure to disk, for later use.
  * Argument file name will get a ".tree" file extension, and in MPI an
  * additional "-[rank]". */
template<int D>
bool FunctionTree_S<D>::saveTree(const string &file) {
    NOT_IMPLEMENTED_ABORT;
//    stringstream fname;
//    fname << file;
//    if (this->isScattered()) {
//        fname << "-" << this->getRankId();
//    }
//    fname << ".tree";
//    ofstream ofs(fname.str().c_str(), ios_base::binary);
//    if (ofs == 0) {
//        MSG_FATAL("Could not open file for writing: " << file);
//    }
//    boost::archive::binary_oarchive oa(ofs);
//    this->purgeGenNodes();
//    oa << *this;
//    return true;
}

/** Read a previously stored tree structure from disk.
  * Argument file name will get a ".tree" file extension, and in MPI an
  * additional "-[rank]". */
template<int D>
bool FunctionTree_S<D>::loadTree(const string &file) {
    NOT_IMPLEMENTED_ABORT;
//    stringstream fname;
//    fname << file;
//    if (node_group.size() > 1 and this->isBuildDistributed()) {
//        fname << "-" << this->getRankId();
//    }
//    fname << ".tree";
//    ifstream ifs(fname.str().c_str(), ios_base::binary);
//    if (not ifs) {
//        return false;
//    }
//    boost::archive::binary_iarchive ia(ifs);
//    ia >> *this;
//    return true;
}

template<int D>
double FunctionTree_S<D>::integrate() const {
    NOT_IMPLEMENTED_ABORT;
    /*    double result = 0.0;
    for (int i = 0; i < this->getRootBox().size(); i++) {
        const FunctionNode<D> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }
#ifdef HAVE_MPI
    result = mpi::all_reduce(node_group, result, std::plus<double>());
#endif
return result;*/
}

template<int D>
double FunctionTree_S<D>::dot(const FunctionTree_S<D> &ket) {
    NOT_IMPLEMENTED_ABORT;
    /*    const FunctionTree_S<D> &bra = *this;
    if (bra.getMRA() != ket.getMRA()) MSG_FATAL("Trees not compatible");
#ifdef HAVE_MPI
    NOT_IMPLEMENTED_ABORT;
//    if (this->isScattered() or rhs.isScattered()) {
//        set<MWNode<D> *> missing;
//        rhs.findMissingInnerProd(*this, missing);
//        rhs.syncNodes(missing);
//    }
#endif
    MRNodeVector nodeTable;
    HilbertIterator<D> it(this);
    it.setReturnGenNodes(false);
    while(it.next()) {
        MRNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const FunctionNode<D> &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MRNode<D> *mrNode = ket.findNode(braNode.getNodeIndex());
        if (mrNode == 0) continue;

        const FunctionNode<D> &ketNode = static_cast<const FunctionNode<D> &>(*mrNode);
        if (braNode.isRootNode()) {
            locResult += braNode.dotScaling(ketNode);
        }
        locResult += braNode.dotWavelet(ketNode);
    }
//#pragma omp critical
    result += locResult;
//    }
#ifdef HAVE_MPI
    NOT_IMPLEMENTED_ABORT;
//    if (this->isScattered()) {
//        return mpi::all_reduce(node_group, result, std::plus<double>());
//    }
#endif
//    this->purgeGenNodes();
//    rhs.purgeGenNodes();
return result;*/
}

template<int D>
double FunctionTree_S<D>::evalf(const double *r) {
    NOT_IMPLEMENTED_ABORT;
/*     double result;
#ifdef HAVE_MPI
    NOT_IMPLEMENTED_ABORT;
//    bool iAmMaster = false;
//    if (this->getRankId() == 0) iAmMaster = true;
//    MRNode<D> &mw_node = this->getNodeOrEndNode(r, 1);
//    FunctionNode<D> &f_node = static_cast<FunctionNode<D> &>(mw_node);
//    if (not f_node.isForeign() or (f_node.isCommon() and iAmMaster)) {
//        result = f_node.asFuncNode().evalf(r);
//        if (not iAmMaster) {
//            node_group.send(0, 0, result);
//        }
//    } else if (iAmMaster) {
//        node_group.recv(mpi::any_source, 0, result);
//    } else {
//        result = 0.0;
//    }
#else
       MRNode<D> &mr_node = this->getNodeOrEndNode(r);
    FunctionNode<D> &f_node = static_cast<FunctionNode<D> &>(mr_node);
    result = f_node.evalf(r);
#endif
    this->deleteGenerated();
    return result;*/
}

template<int D>
void FunctionTree_S<D>::square() {
    NOT_IMPLEMENTED_ABORT
    // Doesn't work in MPI
//    this->purgeGenNodes();
//    for (int i = 0; i < this->endNodeTable.size(); i++) {
//        MWNode<D> &node = *this->endNodeTable[i];
//        node.mwTransform(Reconstruction);
//        node.cvTransform(MWNode<D>::Forward);
//        node.getCoefs() = node.getCoefs().array().square();
//        node.cvTransform(MWNode<D>::Backward);
//        node.mwTransform(Compression);
//    }
//    this->mwTransformUp();
//    this->cropTree();
//    this->squareNorm = this->calcSquareNorm();
}

template<int D>
void FunctionTree_S<D>::power(double d) {
    NOT_IMPLEMENTED_ABORT;
//    this->purgeGenNodes();
//    for (int i = 0; i < this->endNodeTable.size(); i++) {
//        MWNode<D> &node = *this->endNodeTable[i];
//        node.mwTransform(Reconstruction);
//        node.cvTransform(MWNode<D>::Forward);
//        node.getCoefs() = node.getCoefs().array().pow(d);
//        node.cvTransform(MWNode<D>::Backward);
//        node.mwTransform(Compression);
//    }
//    this->mwTransformUp();
//    this->cropTree();
//    this->squareNorm = this->calcSquareNorm();
}

template<int D>
void FunctionTree_S<D>::normalize() {
    NOT_IMPLEMENTED_ABORT;
//    double norm = sqrt(this->getSquareNorm());
//    *this *= (1.0/norm);
}

template<int D>
void FunctionTree_S<D>::orthogonalize(const FunctionTree_S<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
//    this->purgeGenNodes();
//    this->purgeForeignNodes();
//    tree.purgeGenNodes();
//    tree.purgeForeignNodes();
//    double innerProd = this->dot(tree);
//    double norm = tree.getSquareNorm();
//    this->purgeGenNodes();
//    this->purgeForeignNodes();
//    tree.purgeGenNodes();
//    tree.purgeForeignNodes();
//    *this -= (innerProd/norm) * tree;
}

template<int D>
void FunctionTree_S<D>::map(const RepresentableFunction<1> &func) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree_S<D>& FunctionTree_S<D>::operator *=(double c) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree_S<D>& FunctionTree_S<D>::operator *=(const FunctionTree_S<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree_S<D>& FunctionTree_S<D>::operator +=(const FunctionTree_S<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree_S<D>& FunctionTree_S<D>::operator -=(const FunctionTree_S<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template class FunctionTree_S<1> ;
template class FunctionTree_S<2> ;
template class FunctionTree_S<3> ;
