#include "MRGrid.h"
#include "MathUtils.h"
#include "TelePrompter.h"
#include "GridNodeBox.h"
#include "GridNode.h"

using namespace std;

template<int D>
MRGrid<D>::MRGrid(int _order, const GridNodeBox<D> *box) {
    this->order = _order;
    this->kp1 = this->order + 1;
    this->kp1_d = MathUtils::ipow(this->kp1, D);

    if (box != 0) {
        this->rootBox = new GridNodeBox<D>(*box);
    } else {
	NOT_IMPLEMENTED_ABORT
    }

    this->nodesAtDepth.push_back(0);
    initializeRootNodes();
    this->resetEndNodeTable();

    this->maxDepth = 25;
    this->maxScale = this->getRootScale() + this->maxDepth - 1;
}

template<int D>
void MRGrid<D>::initializeRootNodes() {
    GridNode<D> **nodes = this->getRootBox().getNodes();
    assert(nodes != 0);
    for (int i = 0; i < this->getRootBox().getNBoxes(); i++) {
	NodeIndex<D> idx = this->getRootBox().getNodeIndex(i);
        nodes[i] = new GridNode<D>(this, idx);
    }
}

template<int D>
void MRGrid<D>::resetEndNodeTable() {
    NOT_IMPLEMENTED_ABORT
}
    
template<int D>
MRGrid<D>::~MRGrid() {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
const double* MRGrid<D>::getLowerBounds() const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
const double* MRGrid<D>::getUpperBounds() const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int MRGrid<D>::getNNodes(int depth) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int MRGrid<D>::getNLeafNodes(int depth) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int MRGrid<D>::getNQuadPoints(int depth) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int MRGrid<D>::getNQuadPointsPerNode() const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::getQuadraturePoints(Eigen::MatrixXd &roots, const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::getQuadratureWeights(Eigen::VectorXd &weights, const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::getQuadraturePoints(Eigen::MatrixXd &roots) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::getQuadratureWeights(Eigen::VectorXd &weights) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::saveGrid(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::loadGrid(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template class MRGrid<1>;
template class MRGrid<2>;
template class MRGrid<3>;
