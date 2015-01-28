#include "MWProjector.h"
#include "MWNode.h"
#include "GridAdaptor.h"
#include "MRGrid.h"
#include "FunctionTree.h"

using namespace std;

template<int D>
MWProjector<D>::MWProjector() {
    this->tree = 0;
    this->adaptor = 0;
}

template<int D>
MWProjector<D>::MWProjector(GridAdaptor<D> &a) {
    this->tree = 0;
    this->adaptor = &a;
}

template<int D>
MWProjector<D>::~MWProjector() {
    this->adaptor = 0;
    if (this->tree != 0) {
        MSG_ERROR("Projector not properly cleared");
    }
}

template<int D>
void MWProjector<D>::buildTree() {
    println(1, "  == Building tree");
    MRNodeVector nodeTable;
    this->tree->copyEndNodeTable(nodeTable);
    this->tree->clearEndNodeTable();

    int iteration = 1;
    while (nodeTable.size() > 0) {
        calcNodeTable(nodeTable);
        this->tree->calcTreeNorm(&nodeTable);
        if (this->adaptor != 0) {
            NOT_IMPLEMENTED_ABORT;
//            nodeTable = this->adaptor->splitNodeTable(nodeTable);
        } else {
            nodeTable.clear();
        }
        iteration++;
    }
    this->tree->resetEndNodeTable();
}

template<int D>
void MWProjector<D>::calcNodeTable(MRNodeVector &nodeTable) {
    int nNodes = nodeTable.size();
    printout(1, "  -- #" << setw(3) << " Calculated   ");
    printout(1, setw(6) << nNodes << " nodes" << endl);
#pragma omp parallel shared(nodeTable) firstprivate(nNodes)
    {
    #pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = static_cast<MWNode<D> &>(*nodeTable[n]);
            calcWaveletCoefs(node);
        }
    }
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;
