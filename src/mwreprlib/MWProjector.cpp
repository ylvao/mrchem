#include "MWProjector.h"
#include "MWNode.h"
#include "MWAdaptor.h"
#include "MWTree.h"

using namespace std;

template<int D>
MWProjector<D>::MWProjector() {
    this->outTree = 0;
    this->adaptor = 0;
}

template<int D>
MWProjector<D>::MWProjector(MWAdaptor<D> &a) {
    this->outTree = 0;
    this->adaptor = &a;
}

template<int D>
MWProjector<D>::~MWProjector() {
    this->adaptor = 0;
    if (this->outTree != 0) {
        MSG_ERROR("Projector not properly cleared");
    }
}

template<int D>
void MWProjector<D>::buildTree() {
    println(10, " == Building tree");

    MRNodeVector *workVec = this->outTree->copyEndNodeTable();
    MRNodeVector *endVec = this->outTree->getEndNodeTable();
    endVec->clear();

    int iter = 1;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << setw(3) << iter << ": Calculated   ");
        workVec = clearForeignNodes(workVec);
        calcNodeVector(*workVec);
        double norm = sqrt(this->outTree->calcTreeNorm(workVec));
        if (this->adaptor != 0) {
            MRNodeVector splitVec;
            this->adaptor->splitNodeVector(norm, *workVec, splitVec, *endVec);
            NodeIndexSet *splitSet = getNodeIndexSet(splitVec);
            broadcast_index_list<D>(*splitSet);
            this->outTree->splitNodes(*splitSet, workVec);
            delete splitSet;
        } else {
            workVec->clear();
        }
        iter++;
    }
    delete workVec;
    this->outTree->resetEndNodeTable();
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
        calcWaveletCoefs(node);
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
