/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include "GridGenerator.h"
#include "MRGrid.h"
#include "MRNode.h"
#include "TelePrompter.h"

using namespace std;

template<int D>
GridGenerator<D>::GridGenerator(int uScale) : uniformScale(uScale){
    this->grid = 0;
    this->quadOrder = 0;
    this->uniformDepth = -1;
}

template<int D>
GridGenerator<D>::~GridGenerator() {
    if (this->grid != 0) {
        MSG_ERROR("Grid pointer not released");
    }
}

template<int D>
void GridGenerator<D>::generateGrid(MRGrid<D> &outGrid) {
    initGrid(outGrid);
    buildGrid();
    clearGrid();
}

template<int D>
void GridGenerator<D>::initGrid(MRGrid<D> &outGrid) {
    this->grid = &outGrid;
    this->quadOrder = this->grid->getKp1();
    int rootScale = this->grid->getRootScale();
    this->uniformDepth = this->uniformScale - rootScale;
}

template<int D>
void GridGenerator<D>::clearGrid() {
    this->grid = 0;
    this->quadOrder = 0;
    this->uniformDepth = -1;
}

template<int D>
void GridGenerator<D>::buildGrid() {
    println(10, " == Building grid");
    MRNodeVector *nodeVector = this->grid->getEndNodeTable();

    int iter = 1;
    while (nodeVector->size() > 0) {
        int nNodes = nodeVector->size();
        nodeVector = splitNodeVector(nodeVector);
        printout(10, "  -- #" << setw(3) << iter << ": Generated    ");
        printout(10, setw(6) << nNodes << " nodes" << endl);
        iter++;
    }
    this->grid->resetEndNodeTable();
}

template<int D>
MRNodeVector* GridGenerator<D>::splitNodeVector(MRNodeVector *nVec) {
    NodeIndexSet idxSet;
    NodeIndexSet tmpIdx;

    int nNodes = nVec->size();
    for (int n = 0; n < nNodes; n++) {
        MRNode<D> *node = (*nVec)[n];
        if (splitCheck(node)) {
            const NodeIndex<D> *idx = &node->getNodeIndex();
            tmpIdx.insert(idx);
        }
    }
    idxSet.insert(tmpIdx.begin(), tmpIdx.end());
    nVec->clear();
    this->grid->splitNodes(idxSet, nVec);
    return nVec;
}

template<int D>
bool GridGenerator<D>::splitCheck(const MRNode<D> *node) {
    if (node == 0) {
        return false;
    }
    if (node->getDepth() < this->uniformDepth) {
        return true;
    }
    return false;
}

template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
