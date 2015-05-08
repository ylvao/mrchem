#include "MWProjector.h"
#include "MWAdaptor.h"
#include "MWTree.h"
#include "MWNode.h"

using namespace std;

template<int D>
void MWProjector<D>::buildTree(MWTree<D> &outTree) {
    println(10, " == Building tree");

    NodeIndexSet *splitSet = 0;
    MRNodeVector *splitVec = 0;
    MRNodeVector *workVec = outTree.copyEndNodeTable();
    MRNodeVector *endVec = outTree.getEndNodeTable();
    endVec->clear();

    int iter = 0;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << setw(3) << iter << ": Calculated   ");
        workVec = clearForeignNodes(workVec);
        calcNodeVector(*workVec);
        double norm = outTree.calcSquareNorm(workVec);
        if (maxIterReached(iter)) break;
        splitVec = this->adaptor.splitNodeVector(*workVec, endVec);
        splitSet = getNodeIndexSet(*splitVec);
        broadcast_index_list<D>(*splitSet);
        workVec->clear();
        outTree.splitNodes(*splitSet, workVec);
        delete splitSet;
        delete splitVec;
        iter++;
    }
    delete workVec;
    outTree.resetEndNodeTable();
    outTree.calcSquareNorm();
}

template<int D>
void MWProjector<D>::calcNodeVector(MRNodeVector &nodeVec) {
    int nNodes = nodeVec.size();
    printout(10, setw(6) << nNodes << " nodes" << endl);
#pragma omp parallel shared(nodeVec) firstprivate(nNodes)
{
#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        MWNode<D> &node = static_cast<MWNode<D> &>(*nodeVec[n]);
        calcNode(node);
    }
}
}

template<int D>
MRNodeVector* MWProjector<D>::clearForeignNodes(MRNodeVector *oldVec) const {
    MRNodeVector *newVec = new MRNodeVector;
    for (int i = 0; i < oldVec->size(); i++) {
        MRNode<D> *node = (*oldVec)[i];
        if (node == 0) {
            continue;
        }
        if (not node->isForeign()) {
            newVec->push_back(node);
        }
    }
    delete oldVec;
    return newVec;
}

template<int D>
NodeIndexSet* MWProjector<D>::getNodeIndexSet(const MRNodeVector &nodeVec) const {
    NodeIndexSet *idxSet = new NodeIndexSet;
    for (int i = 0; i < nodeVec.size(); i++) {
        const NodeIndex<D> &idx = nodeVec[i]->getNodeIndex();
        idxSet->insert(&idx);
    }
    return idxSet;
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;
