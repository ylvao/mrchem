/**
 *  \date Oct 12, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 */

#include "FunctionTree_S.h"
#include "MWNode.h"
#include "MWTree.h"
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
          lastNode(0),
          mwTree_p(0),
          tree_S_array(0),
          maxNodes(maxNumberOfNodes),
          sizeTreeMeta(0),
          sizeNode(0) {
    //The first part of the Tree is filled with metadata; reserved size:
    this->sizeTreeMeta = (sizeof(FunctionTree<D>)+7)/sizeof(double);
    //The dynamical part of the tree is filled with nodes of size:
    this->sizeNode = (sizeof(ProjectedNode<D>)+7)/sizeof(double);
    cout<<"SizeNode "<<sizeof(ProjectedNode<D>)<<endl;

    //Tree is defined as array of doubles, because C++ does not like void malloc
    this->tree_S_array = new double[this->sizeTreeMeta + maxNumberOfNodes*this->sizeNode];
    this->lastNode = (ProjectedNode<D>*) (this->tree_S_array + this->sizeTreeMeta);

    this->mwTree_p = new (this->tree_S_array) MWTree<D>(mra);//put a MWTree at start of array
    //this->funcTree_p = new FunctionTree<D>(mra);//put a MWTree at start of array

    ProjectedNode<D> *node_p = allocNodes(this->mwTree_p->getRootBox().size());

    NodeBox<D> &rBox = this->mwTree_p->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        const NodeIndex<D> &nIdx = rBox.getNodeIndex(rIdx);
        roots[rIdx] = new ProjectedNode<D>(getTree(), nIdx);
        node_p++;
    }
    this->mwTree_p->resetEndNodeTable();
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
        cout << "new size " << this->nNodes << endl;
        return this->lastNode - nAlloc*this->sizeNode;
    }
}
/** FunctionTree_S destructor. */
template<int D>
FunctionTree_S<D>::~FunctionTree_S() {
    MWNode<D> **roots = this->mwTree_p->getRootBox().getNodes();
    for (int i = 0; i < this->mwTree_p->getRootBox().size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
        node->~ProjectedNode();
        roots[i] = 0;
    }
    this->mwTree_p->~MWTree();
    delete[] this->tree_S_array;
}

template class FunctionTree_S<1> ;
template class FunctionTree_S<2> ;
template class FunctionTree_S<3> ;
