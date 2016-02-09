/**
 *  \date Oct 12, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 */

#include "FunctionTree_S.h"
#include "MWNode.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
template<int D>
FunctionTree_S<D>::FunctionTree_S(const MultiResolutionAnalysis<D> &mra, int maxNumberOfNodes)
        : nNodes(0),
          maxNodes(maxNumberOfNodes) {
    //The first part of the Tree is filled with metadata; reserved size:
    sizeTreeMeta = (sizeof(FunctionTree<D>)+7)/sizeof(double);
    //The dynamical part of the tree is filled with nodes of size:
    sizeNode = (sizeof(ProjectedNode<D>)+7)/sizeof(double);
    cout<<"SizeNode "<<sizeof(ProjectedNode<D>)<<endl;

    //Tree is defined as array of doubles, because C++ does not like void malloc
    this->tree_S_array = new double[this->sizeTreeMeta + maxNumberOfNodes*this->sizeNode];

    this->lastNode = (ProjectedNode<D>*) (this->tree_S_array + this->sizeTreeMeta + this->nNodes*this->sizeNode);

    this->funcTree_p = new (this->tree_S_array) FunctionTree<D>(mra);//put a MWTree at start of array
//    this->funcTree_p = new FunctionTree<D>(mra);//put a MWTree at start of array

    this->lastNode = allocNodes(this->funcTree_p->getRootBox().size());

    for (int rIdx = 0; rIdx < this->funcTree_p->getRootBox().size(); rIdx++) {
        const NodeIndex<D> &nIdx = this->funcTree_p->getRootBox().getNodeIndex(rIdx);
        MWNode<D> *fNode  = new (this->lastNode) ProjectedNode<D>(*((FunctionTree<D>*) this->funcTree_p), nIdx);
//        MWNode<D> *fNode  = new ProjectedNode<D>(*((FunctionTree<D>*) this->funcTree_p), nIdx);
        this->funcTree_p->getRootBox().setNode(rIdx, &fNode);
        this->lastNode++;
    }
    this->funcTree_p->resetEndNodeTable();
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* FunctionTree_S<D>::allocNodes(int nAlloc) {
    this->nNodes += nAlloc;
    if (this->nNodes > this->maxNodes){
        this->nNodes -= nAlloc;
        return 0;
    } else {
        this->lastNode += nAlloc*this->sizeNode;
        cout<<"new size "<<this->nNodes<<endl;
        return this->lastNode-nAlloc*this->sizeNode;
    }
}
/** FunctionTree_S destructor. */
template<int D>
FunctionTree_S<D>::~FunctionTree_S() {
    delete[] this->tree_S_array;
    //delete this->funcTree_p;
}

template class FunctionTree_S<1> ;
template class FunctionTree_S<2> ;
template class FunctionTree_S<3> ;
