#include "MWProjector.h"
#include "MWNode.h"
#include "GridAdaptor.h"
#include "MRGrid.h"
#include "FunctionTree.h"

using namespace std;

template<int D>
MWProjector<D>::MWProjector() {
    this->outTree = 0;
    this->adaptor = 0;
}

template<int D>
MWProjector<D>::MWProjector(GridAdaptor<D> &a) {
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
    MRNodeVector nodeTable;
    this->outTree->copyEndNodeTable(nodeTable);
    this->outTree->clearEndNodeTable();

    int iteration = 1;
    while (nodeTable.size() > 0) {
        calcNodeTable(nodeTable);
        this->outTree->calcTreeNorm(&nodeTable);
        if (this->adaptor != 0) {
            NOT_IMPLEMENTED_ABORT;
//            nodeTable = this->adaptor->splitNodeTable(nodeTable);
        } else {
            nodeTable.clear();
        }
        iteration++;
    }
    this->outTree->resetEndNodeTable();
}

template<int D>
void MWProjector<D>::calcNodeTable(MRNodeVector &nodeTable) {
    int nNodes = nodeTable.size();
    printout(10, "  -- #  1: Calculated   ");
    printout(10, setw(6) << nNodes << " nodes" << endl);
#pragma omp parallel shared(nodeTable) firstprivate(nNodes)
    {
    #pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = static_cast<MWNode<D> &>(*nodeTable[n]);
            if (not node.isForeign()) {
                calcWaveletCoefs(node);
            }
        }
    }
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;
