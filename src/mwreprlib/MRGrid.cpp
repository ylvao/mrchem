#include "MRGrid.h"
#include "MathUtils.h"
#include "TelePrompter.h"
#include "NodeBox.h"
#include "GridNode.h"
#include "LebesgueIterator.h"

using namespace Eigen;
using namespace std;

template<int D>
MRGrid<D>::MRGrid(const BoundingBox<D> *box, int k) : MRTree<D>(box, k) {
    initializeRootNodes();
    this->resetEndNodeTable();
}

template<int D>
void MRGrid<D>::initializeRootNodes() {
    for (int i = 0; i < this->getNRootNodes(); i++) {
        NodeIndex<D> nIdx = this->getRootBox().getNodeIndex(i);
        MRNode<D> *root = new GridNode<D>(*this, nIdx);
        this->rootBox->setNode(i, &root);
    }
}

template<int D>
MRGrid<D>::~MRGrid() {
}

template<int D>
void MRGrid<D>::clear() {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int MRGrid<D>::countQuadPoints(int depth) {
    int nNodes = this->countLeafNodes(depth);
    int ptsPerNode = this->tDim*this->kp1_d;
    return nNodes*ptsPerNode;
}

template<int D>
void MRGrid<D>::getQuadPoints(Eigen::MatrixXd &gridPoints) {
    int nPoints = this->getKp1_d();
    int nNodes = this->endNodeTable.size();
    int totPoints = nNodes*nPoints;

    gridPoints = MatrixXd::Zero(totPoints,D);
    MatrixXd nodePoints = MatrixXd(nPoints,D);
    for (int n = 0; n < nNodes; n++) {
        GridNode<D> &node = static_cast<GridNode<D> &>(*this->endNodeTable[n]);
        node.getExpandedPoints(nodePoints);
        gridPoints.block(n*nPoints, 0, nPoints, D) = nodePoints;
    }
}

template<int D>
void MRGrid<D>::getQuadWeights(Eigen::VectorXd &gridWeights) {
    int nWeights = this->getKp1_d();
    int nNodes = this->endNodeTable.size();
    int totWeights = nNodes*nWeights;

    gridWeights = VectorXd::Zero(totWeights);
    VectorXd nodeWeights = VectorXd(nWeights);

    for (int n = 0; n < nNodes; n++) {
        GridNode<D> &node = static_cast<GridNode<D> &>(*this->endNodeTable[n]);
        node.getExpandedWeights(nodeWeights);
        for (int i = 0; i < nWeights; i++) {
            gridWeights(n*nWeights + i) = nodeWeights(i);
        }
    }
}

template<int D>
bool MRGrid<D>::saveTree(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
bool MRGrid<D>::loadTree(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template class MRGrid<1>;
template class MRGrid<2>;
template class MRGrid<3>;
