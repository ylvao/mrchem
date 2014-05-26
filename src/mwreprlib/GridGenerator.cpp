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
#include "GridNode.h"
#include "TelePrompter.h"

using namespace std;

template<int D>
GridGenerator<D>::~GridGenerator() {
    if (this->grid != 0) {
        THROW_ERROR("Grid pointer not released");
    }
}

template<int D>
void GridGenerator<D>::generateGrid(MRGrid<D> &outGrid) { 
    println(0, " == Generating grid ");
    this->grid = &outGrid;
    GridNodeVector nodeTable;
    this->grid->copyEndNodeTable(nodeTable);
    this->grid->clearEndNodeTable();

    int iteration = 1;
    while (nodeTable.size() > 0) {
        int nNodes = nodeTable.size();
        splitNodeTable(nodeTable);
        println(0, "  -- #" << setw(3) << iteration << ": Generated    "
                << setw(6) << nNodes << " nodes");
        iteration++;
    }
    this->grid->resetEndNodeTable();
    this->grid = 0;
}

template<int D>
void GridGenerator<D>::splitNodeTable(GridNodeVector &nodeTable) { 
    NodeIndexSet idxSet;
    NodeIndexSet tmpIdx;

    // nodeTable contains ALL endNodes, including foreign
    int nNodes = nodeTable.size();
#pragma omp parallel firstprivate(nNodes) \
    private(tmpIdx, tmpEndNodes)
    {
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            GridNode<D> *node = nodeTable[n];
            if (splitCheck(node)) {
                const NodeIndex<D> *idx = &node->getNodeIndex();
                tmpIdx.insert(idx);
            }
        }
#pragma omp critical
        {
            idxSet.insert(tmpIdx.begin(), tmpIdx.end());
        }
    }
    nodeTable.clear();
    this->grid->yieldChildren(nodeTable, idxSet);
}

template<int D>
bool GridGenerator<D>::splitCheck(GridNode<D> *node) {
    if (node == 0) {
	return false;
    }
    if (node->getDepth() < this->uniform) {
	return true;
    }
    return false;
}


template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
