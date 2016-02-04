/**
 *  \date April 20, 2010
 *  CTCC, University of Troms
 *
 */

#include "MWTree.h"
#include "MultiResolutionAnalysis.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

/** MWTree constructor.
  * Creates an empty tree object. Node construction and assignment of most of
  * the parameters are done in derived classes. */
template<int D>
MWTree<D>::MWTree(const MultiResolutionAnalysis<D> &mra)
        : MRTree<D>(mra.getWorldBox()),
          MRA(mra),
          order(mra.getOrder()),
          kp1_d(MathUtils::ipow(mra.getOrder() + 1, D)),
          squareNorm(-1.0) {
    allocWorkMemory();
}

/** MWTree copy constructor.
  * Takes the parameters of the input tree, not it's data */
/*
template<int D>
MWTree<D>::MWTree(const MRGrid<D> &grid, int type) : MRTree<D>(grid) {
    this->squareNorm = 0.0;
    this->scalingType = type;

    setupScalingBasis(this->scalingType, this->order);
    setupFilters(this->scalingType, this->order);

    allocWorkMemory();
}
*/

/** MWTree copy constructor.
  * Takes the parameters of the input tree, not it's data */
template<int D>
MWTree<D>::MWTree(const MWTree<D> &tree)
        : MRTree<D>(tree),
          MRA(tree.MRA),
          order(tree.order),
          kp1_d(tree.kp1_d),
          squareNorm(-1.0) {
    allocWorkMemory();
}

/** MWTree destructor. */
template<int D>
MWTree<D>::~MWTree() {
    freeWorkMemory();
}

/** Allocate work memory of the tree, for use in mwTransform and the like. */
template<int D>
void MWTree<D>::allocWorkMemory() {
    this->tmpCoefs = new Eigen::MatrixXd *[this->nThreads];
    this->tmpVector = new Eigen::VectorXd *[this->nThreads];
    this->tmpMWCoefs = new Eigen::VectorXd *[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->tmpCoefs[i] = new Eigen::MatrixXd(this->order + 1, D);
        this->tmpVector[i] = new Eigen::VectorXd(this->kp1_d);
        this->tmpMWCoefs[i] = new Eigen::VectorXd(this->kp1_d * (1 << D));
    }
}

/** Deallocate work memory */
template<int D>
void MWTree<D>::freeWorkMemory() {
    for (int i = 0; i < this->nThreads; i++) {
        delete this->tmpCoefs[i];
        delete this->tmpVector[i];
        delete this->tmpMWCoefs[i];
    }
    delete[] this->tmpCoefs;
    delete[] this->tmpVector;
    delete[] this->tmpMWCoefs;
}


template<int D>
double MWTree<D>::estimateError(bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
//    double error = 0.0;
//    for (int i = 0; i < this->getNEndNodes(); i++) {
//        MWNode<D> &node = getEndMWNode(i);
//        error += node.estimateError(absPrec);
//    }
//#ifdef HAVE_MPI
//    error = mpi::all_reduce(node_group, error, std::plus<double>());
//#endif
//    return error;
}

/** Calculate the squared norm of a function represented as a tree.
  *
  * Norm is calculated using endNodes only, but if your endNodeTable is
  * incomplete (e.g. within growTree), the missing nodes must be given in the
  * input work vector. Involves an MPI reduction operation. */
template<int D>
void MWTree<D>::calcSquareNorm(const MRNodeVector *nodeVec)  {
    double treeNorm = 0.0;

    int nNodes = 0;
    if (nodeVec != 0) nNodes = nodeVec->size();

    bool doCommon = false;
    if (this->getRankId() == 0) doCommon = true;

    for (int n = 0; n < nNodes; n++) {
        const MWNode<D> &node = static_cast<const MWNode<D> &>(*(*nodeVec)[n]);
        if (node.isLocal() or (node.isCommon() and doCommon)) {
            assert(node.hasCoefs());
            treeNorm += node.getSquareNorm();
        }
    }
    for (int n = 0; n < this->getNEndNodes(); n++) {
        const MWNode<D> &node = getEndMWNode(n);
        if (node.isLocal() or (node.isCommon() and doCommon)) {
            assert(node.hasCoefs());
            treeNorm += node.getSquareNorm();
        }
    }
#ifdef HAVE_MPI
    this->squareNorm = mpi::all_reduce(node_group, treeNorm, std::plus<double>());
#else
    this->squareNorm = treeNorm;
#endif
}

/** Reduce the accuracy of the tree by deleting nodes
  * which have a higher precision than the requested precison.
  * By default, the relative precision of the tree is used. */
template<int D>
void MWTree<D>::crop(double thrs, bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
}

//template<int D>
//void MWTree<D>::cropTree(double prec, bool absPrec) {
//    set<const NodeIndex<D> *, NodeIndexComp<D> > cropNodes;
//    for (int i = 0; i < this->rootBox.getNBoxes(); i++) {
//        MWNode<D> &rootNode = getRootMWNode(i);
//        if (this->isScattered()) {
//            rootNode.cropChildren(prec, &cropNodes);
//        } else {
//            rootNode.cropChildren(prec);
//        }
//    }
//    if (this->isScattered()) {
//        broadcastIndexList(cropNodes);
//        typename set<const NodeIndex<D> *,
//                NodeIndexComp<D> >::reverse_iterator it;
//        for (it = cropNodes.rbegin(); it != cropNodes.rend(); it++) {
//            MWNode<D> *node = this->findNode(**it);
//            if (node != 0) {
//                node->deleteChildren();
//            }
//        }
//    }
//    resetEndNodeTable();
//    this->squareNorm = calcSquareNorm();
//}

template<int D>
void MWTree<D>::mwTransform(int type, bool overwrite) {
    switch(type) {
    case TopDown:
        mwTransformDown(overwrite);
        break;
    case BottomUp:
        mwTransformUp(overwrite);
        break;
    default:
        MSG_FATAL("Invalid wavelet transform");
    }
}

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
  * thus purifying all coefficients. Option to overwrite or add up existing
  * coefficients of BranchNodes (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformUp(bool overwrite) {
    vector<MRNodeVector > nodeTable;
    this->makeLocalNodeTable(nodeTable, true);

    int start = nodeTable.size() - 2;
    for (int n = start; n >= 0; n--) {
//        set<MRNode<D> *> missing;
//        this->findMissingChildren(nodeTable[n], missing);
        //communicate missing
//        this->syncNodes(missing);
        int nNodes = nodeTable[n].size();
#pragma omp parallel firstprivate(nNodes, overwrite) shared(nodeTable)
        {
#pragma omp for schedule(guided)
            for (int i = 0; i < nNodes; i++) {
                MWNode<D> &node = static_cast<MWNode<D> &>(*nodeTable[n][i]);
                if ((not node.isForeign()) and node.isBranchNode()) {
                    node.reCompress(overwrite);
                    node.calcNorms();
                }
            }
        }
    }
}

/** Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs
  * on coarser scales, starting at the rootNodes. Option to overwrite or add up
  * existing scaling coefficients (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformDown(bool overwrite) {
    NOT_IMPLEMENTED_ABORT;
//    this->purgeForeignNodes();
//    if (isScattered()) {
//        vector<MWNodeVector > nodeTable;
//        this->makeLocalNodeTable(nodeTable);

//        vector<MWNodeVector > parentTable;
//        this->makeNodeTable(parentTable);

//        set<MWNode<D> *> missing;
//        for (int i = 1; i < nodeTable.size(); i++) {
//            missing.clear();
//            findMissingParents(nodeTable[i], missing);
//            //communicate missing
//            syncNodes(missing);
//            for (unsigned int j = 0; j < parentTable[i-1].size(); j++) {
//                MWNode<D> &parent = *parentTable[i-1][j];
//                if (not parent.hasCoefs() or parent.isLeafNode()) {
//                    continue;
//                }
//                parent.giveChildrenScaling(overwrite);
//                for (int k = 0; k < parent.getNChildren(); k++) {
//                    MWNode<D> &child = parent.getMWChild(k);
//                    if (child.isForeign()) {
//                        child.clearCoefs();
//                        child.clearNorms();
//                    }
//                }
//            }
//        }
//    } else {
//        vector<MWNodeVector > nodeTable;
//        makeNodeTable(nodeTable);
//#pragma omp parallel shared(nodeTable)
//        {
//            for (int n = 0; n < nodeTable.size(); n++) {
//                int n_nodes = nodeTable[n].size();
//#pragma omp for schedule(guided)
//                for (int i = 0; i < n_nodes; i++) {
//                    MWNode<D> &node = *nodeTable[n][i];
//                    if (node.isBranchNode()) {
//                        node.giveChildrenScaling(overwrite);
//                    }
//                }
//            }
//        }
//        this->purgeForeignNodes();
//    }
}

/** Traverse tree and set all nodes to zero.
  *
  * Keeps the node structure of the tree, even though the zero function is
  * representable at depth zero. Use cropTree to remove unnecessary nodes.*/
template<int D>
void MWTree<D>::setZero() {
    NOT_IMPLEMENTED_ABORT;
//    HilbertIterator<D> it(this);
//    while(it.next()) {
//        MWNode<D> &node = static_cast<MWNode<D> &>(it.getNode());
//        if (not node.isForeign()) {
//            node.getCoefs().setZero();
//            node.calcNorms();
//        }
//    }
//    this->squareNorm = 0.0;
}

template class MWTree<1>;
template class MWTree<2>;
template class MWTree<3>;
