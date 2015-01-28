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
    println(1, " == Building grid");
    MRNodeVector nodeTable;
    this->grid->copyEndNodeTable(nodeTable);
    this->grid->clearEndNodeTable();

    int iteration = 1;
    while (nodeTable.size() > 0) {
        int nNodes = nodeTable.size();
        splitNodeTable(nodeTable);
        println(1, "  -- #" << setw(3) << iteration << ": Generated    "
                << setw(6) << nNodes << " nodes");
        iteration++;
    }
    this->grid->resetEndNodeTable();
}

template<int D>
void GridGenerator<D>::splitNodeTable(MRNodeVector &nodeTable) {
    NodeIndexSet idxSet;
    NodeIndexSet tmpIdx;

    int nNodes = nodeTable.size();
    for (int n = 0; n < nNodes; n++) {
        MRNode<D> *node = nodeTable[n];
        if (splitCheck(node)) {
            const NodeIndex<D> *idx = &node->getNodeIndex();
            tmpIdx.insert(idx);
        }
    }
    idxSet.insert(tmpIdx.begin(), tmpIdx.end());
    nodeTable.clear();
    this->grid->yieldChildren(nodeTable, idxSet);
}

template<int D>
bool GridGenerator<D>::splitCheck(const MRNode<D> *node) {
    if (node == 0) {
        return false;
    }
    if (node->getDepth() < this->uniformDepth) {
        //println(0, "uniform depth split " << node->getDepth() << " < " << this->uniformDepth);
        return true;
    }
    return false;
}

template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
