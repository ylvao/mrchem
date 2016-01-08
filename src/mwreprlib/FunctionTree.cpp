/**
 *  \date Oct 12, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 */

#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
template<int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra)
        : MWTree<D> (mra) {
    const double *lB = this->rootBox.getLowerBounds();
    const double *uB = this->rootBox.getUpperBounds();
    this->setBounds(lB, uB);
    for (int rIdx = 0; rIdx < this->rootBox.size(); rIdx++) {
        const NodeIndex<D> &nIdx = this->rootBox.getNodeIndex(rIdx);
        MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
        this->rootBox.setNode(rIdx, &root);
    }
    this->resetEndNodeTable();
    this->calcSquareNorm();
}

/** FunctionTree copy constructor.
  * Copy polynomial order and type, as well as the world box from the
  * given tree, but only at root scale. Initializes the function to zero.
  * Use = operator to copy data.*/
template<int D>
FunctionTree<D>::FunctionTree(const MWTree<D> &tree) : MWTree<D> (tree) {
    NOT_IMPLEMENTED_ABORT;
//    const double *lB = this->rootBox->getLowerBounds();
//    const double *uB = this->rootBox->getUpperBounds();
//    this->setBounds(lB, uB);
//    for (int rIdx = 0; rIdx < this->getNRootNodes(); rIdx++) {
//        const NodeIndex<D> &nIdx = this->rootBox->getNodeIndex(rIdx);
//        MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
//        this->rootBox->setNode(rIdx, &root);
//    }
//    this->resetEndNodeTable();
}

/** FunctionTree copy constructor.
  * Copy polynomial order and type, as well as the world box from the
  * given tree, but only at root scale. Initializes the function to zero.
  * Use = operator to copy data.*/
template<int D>
FunctionTree<D>::FunctionTree(const FunctionTree<D> &tree) : MWTree<D> (tree) {
    NOT_IMPLEMENTED_ABORT;
//    const double *lB = this->rootBox->getLowerBounds();
//    const double *uB = this->rootBox->getUpperBounds();
//    this->setBounds(lB, uB);
//    for (int rIdx = 0; rIdx < this->getNRootNodes(); rIdx++) {
//        const NodeIndex<D> &nIdx = this->rootBox->getNodeIndex(rIdx);
//        MRNode<D> *root = new ProjectedNode<D>(*this, nIdx);
//        this->rootBox->setNode(rIdx, &root);
//    }
//    this->resetEndNodeTable();
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

/** FunctionTree destructor. */
template<int D>
FunctionTree<D>::~FunctionTree() {
    MRNode<D> **rootNodes = this->rootBox.getNodes();
    for (int i = 0; i < this->rootBox.size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(rootNodes[i]);
        if (node != 0) delete node;
        rootNodes[i] = 0;
    }
}

/** Leaves the tree inn the same state as after construction*/
template<int D>
void FunctionTree<D>::clear() {
    NOT_IMPLEMENTED_ABORT;
}

/** Loop through endNodeTable and recursively clear all GenNode coefficients.
  * Includes a static cast of endNodes from MWNode to FunctionNode*/
template<int D>
void FunctionTree<D>::clearGenNodes() {
    NOT_IMPLEMENTED_ABORT;
//    for (unsigned int i = 0; i < this->endNodeTable.size(); i++) {
//        FunctionNode<D> *node =
//                static_cast<FunctionNode<D> *>(this->endNodeTable[i]);
//        node->clearGenerated();
//    }
}

/** Loop through endNodeTable and recursively delete all GenNodes.
  * Includes a static cast of endNodes from MWNode to FunctionNode*/
template<int D>
void FunctionTree<D>::purgeGenNodes() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nNodes == 0) {
//        return;
//    }
//    int nEnd = this->endNodeTable.size();
//#pragma omp parallel firstprivate(nEnd)
//    {
//#pragma omp for schedule(guided)
//        for (int n = 0; n < nEnd; n++) {
//            this->getEndNode(n).purgeGenerated();
//        }
//    }
}

/** Write the tree structure to disk, for later use.
  * Argument file name will get a ".tree" file extension, and in MPI an
  * additional "-[rank]". */
template<int D>
bool FunctionTree<D>::saveTree(const string &file) {
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
bool FunctionTree<D>::loadTree(const string &file) {
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
FunctionNode<D>& FunctionTree<D>::getRootFuncNode(int rIdx) {
    return static_cast<FunctionNode<D> &>(this->getRootNode(rIdx));
}

template<int D>
const FunctionNode<D>& FunctionTree<D>::getRootFuncNode(int rIdx) const {
    return static_cast<const FunctionNode<D> &>(this->getRootNode(rIdx));
}

template<int D>
FunctionNode<D>& FunctionTree<D>::getRootFuncNode(const NodeIndex<D> &nIdx) {
    return static_cast<FunctionNode<D> &>(this->getRootNode(nIdx));
}

template<int D>
const FunctionNode<D>& FunctionTree<D>::getRootFuncNode(const NodeIndex<D> &nIdx) const {
    return static_cast<const FunctionNode<D> &>(this->getRootNode(nIdx));
}

template<int D>
double FunctionTree<D>::integrate() {
    NOT_IMPLEMENTED_ABORT;
//    double result = 0.0;
//    for (int i = 0; i < this->getNRootNodes(); i++) {
//        FunctionNode<D> &fNode = getRootFuncNode(i);
//        result += fNode.integrate();
//    }
//#ifdef HAVE_MPI
//    result = mpi::all_reduce(node_group, result, std::plus<double>());
//#endif
//    return result;
}

template<int D>
double FunctionTree<D>::dot(FunctionTree<D> &ket) {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->checkCompatible(rhs)) {
//        MSG_FATAL("Trees not compatible");
//    }
//#ifdef HAVE_MPI
//    if (this->isScattered() or rhs.isScattered()) {
//        set<MWNode<D> *> missing;
//        rhs.findMissingInnerProd(*this, missing);
//        rhs.syncNodes(missing);
//    }
//#endif
//    MWNodeVector nodeTable;
//    HilbertIterator<D> it(this);
//    while(it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        nodeTable.push_back(&node);
//    }
//    int n_nodes = nodeTable.size();
//    double result = 0.0;
//    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
//    for (int n = 0; n < n_nodes; n++) {
//        FunctionNode<D> *nodeA = static_cast<FunctionNode<D> *>
//                (nodeTable[n]);
//        FunctionNode<D> *nodeB = static_cast<FunctionNode<D> *>
//                (rhs.findNode(nodeA->getNodeIndex()));
//        if (nodeB == 0) {
//        continue;
//        }
//            if (nodeA->isRoot()) {
//            locResult += nodeA->scalingInnerProduct(*nodeB);
//        }
//        locResult += nodeA->waveletInnerProduct(*nodeB);
//    }
//#pragma omp critical
//    result += locResult;
//    }
//#ifdef HAVE_MPI
//    if (this->isScattered()) {

//        return mpi::all_reduce(node_group, result, std::plus<double>());
//    }
//#endif
//    this->purgeGenNodes();
//    rhs.purgeGenNodes();
//    return result;
}

template<int D>
double FunctionTree<D>::evalf(const double *r) const {
    NOT_IMPLEMENTED_ABORT;
//    double val;
//    FunctionTree &tree = const_cast<FunctionTree &> (*this);

//    if (this->isScattered()) {
//#ifdef HAVE_MPI

//        MWNode<D> *pnode = tree.findNode(r);
//        if (not pnode->isForeign()) {
//            MWNode<D> &node = pnode->getNode(r);
//            val = node.evalf(r);
//            if (this->getRankId() != 0) {
//                node_group.send(0, 0, val);
//            }
//        } else if (this->getRankId() == 0) {
//            node_group.recv(mpi::any_source, 0, val);
//        } else {
//            val = 0.0;
//        }
//#else
//        MSG_FATAL("Calling evalf() on scattered tree without MPI enabled!");
//#endif
//    } else {
//        MWNode<D> &node = tree.getNode(r);
//        val = node.evalf(r);
//    }
//    return val;
}

template<int D>
void FunctionTree<D>::square() {
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
void FunctionTree<D>::power(double d) {
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
void FunctionTree<D>::normalize() {
    NOT_IMPLEMENTED_ABORT;
//    double norm = sqrt(this->getSquareNorm());
//    *this *= (1.0/norm);
}

template<int D>
void FunctionTree<D>::orthogonalize(const FunctionTree<D> &tree) {
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
void FunctionTree<D>::map(const RepresentableFunction<1> &func) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator *=(double c) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator *=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator +=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator -=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template class FunctionTree<1> ;
template class FunctionTree<2> ;
template class FunctionTree<3> ;
