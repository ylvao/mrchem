#include "MRGrid.h"
#include "MathUtils.h"
#include "TelePrompter.h"
#include "GridNodeBox.h"
#include "GridNode.h"
#include "LebesgueIterator.h"

using namespace Eigen;
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

    this->maxDepth = 25;
    this->maxScale = this->getRootScale() + this->maxDepth - 1;

    this->nodesAtDepth.push_back(0);
    initializeRootNodes();
    this->resetEndNodeTable();
}

template<int D>
void MRGrid<D>::initializeRootNodes() {
    for (int i = 0; i < this->getRootBox().getNBoxes(); i++) {
	NodeIndex<D> idx = this->getRootBox().getNodeIndex(i);
        GridNode<D> *root = new GridNode<D>(this, idx.getScale(), idx.getTranslation());
	this->rootBox->setNode(i, &root);
    }
}

template<int D>
MRGrid<D>::~MRGrid() {
    this->endNodeTable.clear();
    delete this->rootBox;
}

template<int D>
void MRGrid<D>::incrementNodeCount(int scale) {
    int depth = scale - this->getRootScale();
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
        for (int i = 0; i < depth - n; i++) {
            this->nodesAtDepth.push_back(0);
        }
    }
    int nodes = this->nodesAtDepth[depth];
    nodes++;
    this->nodesAtDepth[depth] = nodes;
}

template<int D>
void MRGrid<D>::decrementNodeCount(int scale) {
    unsigned int depth = scale - this->getRootScale();
    if (depth >= this->nodesAtDepth.size()) {
        THROW_ERROR("Depth out of bounds!");
    }
    int nodes = this->nodesAtDepth[depth];
    nodes--;
    if (nodes < 0) {
        THROW_ERROR("Number of nodes cannot be negative.");
    }
    this->nodesAtDepth[depth] = nodes;
    if (nodes == 0) { 
        this->nodesAtDepth.pop_back();
    }
}

template<int D>
int MRGrid<D>::getNodesAtDepth(int depth) const {
    if (depth < 0) {
	THROW_ERROR("Negative depth");
    }
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
	return 0;
    }
    return this->nodesAtDepth[depth];
}

template<int D>
void MRGrid<D>::clearEndNodeTable() {
    this->endNodeTable.clear();
}

template<int D>
void MRGrid<D>::resetEndNodeTable() {
    this->endNodeTable.clear();
    LebesgueIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
        if (node.isEndNode()) {
	    this->endNodeTable.push_back(&node);
        }
    }
}

template<int D>
void MRGrid<D>::copyEndNodeTable(GridNodeVector &outTable) {
    for (int i = 0; i < this->endNodeTable.size(); i++) {
	GridNode<D> *node = this->endNodeTable[i];
	outTable.push_back(node);
    }
}
    
template<int D>
int MRGrid<D>::getNNodes(int depth) const {
    int nScales = this->nodesAtDepth.size();
    if (depth > nScales) {
	return 0;
    }
    if (depth >= 0) {
	return this->nodesAtDepth[depth];
    }
    int totNodes = 0;
    for (int i = 0; i < nScales; i++) {
	totNodes += this->nodesAtDepth[i];
    }
    return totNodes; 
}

template<int D>
int MRGrid<D>::countBranchNodes(int depth) {
    int nNodes = 0;
    LebesgueIterator<D> it(this);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
	if (node.getDepth() == depth or depth < 0) {
            if (node.isBranchNode()) {
		nNodes++;
	    }
	}
    }
    return nNodes;
}

template<int D>
int MRGrid<D>::countLeafNodes(int depth) {
    int nNodes = 0;
    LebesgueIterator<D> it(this);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
	if (node.getDepth() == depth or depth < 0) {
            if (node.isLeafNode()) {
		nNodes++;
	    }
	}
    }
    return nNodes;
}

template<int D>
int MRGrid<D>::countQuadPoints(int depth) {
    int nNodes = countLeafNodes(depth);
    int ptsPerNode = this->tDim*this->kp1_d;
    return nNodes*ptsPerNode;
}

template<int D>
void MRGrid<D>::getQuadPoints(Eigen::MatrixXd &gridPoints) {
    int tDim = getTDim();
    int kp1_d = getKp1_d();
    int nPoints = tDim*kp1_d;
    int nNodes = this->endNodeTable.size();
    int totPoints = nPoints*nNodes;

    gridPoints = MatrixXd::Zero(totPoints,D);
    MatrixXd nodePoints = MatrixXd(nPoints,D);
    for (int n = 0; n < nNodes; n++) {
        GridNode<D> &node = *this->endNodeTable[n];
	node.createChildren();
	for (int cIdx = 0; cIdx < tDim; cIdx++) {
	    GridNode<D> &child = *node.getChild(cIdx);
	    child.getExpandedPoints(nodePoints);
	    int inpos = n*nPoints + cIdx*kp1_d;
	    gridPoints.block(inpos, 0, kp1_d, D) = nodePoints;
	}
	node.deleteChildren();
    }
}

template<int D>
void MRGrid<D>::getQuadWeights(Eigen::VectorXd &gridWeights) {
    int tDim = getTDim();
    int kp1_d = getKp1_d();
    int nWeights = tDim*kp1_d;
    int nNodes = this->endNodeTable.size();
    int totWeights = nWeights*nNodes;

    gridWeights = VectorXd::Zero(totWeights);
    VectorXd nodeWeights = VectorXd(nWeights);

    for (int n = 0; n < nNodes; n++) {
        GridNode<D> &node = *this->endNodeTable[n];
	node.createChildren();
	for (int cIdx = 0; cIdx < tDim; cIdx++) {
	    GridNode<D> &child = *node.getChild(cIdx);
	    child.getExpandedWeights(nodeWeights);
	    int inpos = n*nWeights + cIdx*kp1_d;
	    gridWeights.segment(inpos, kp1_d) = nodeWeights;
	}
	node.deleteChildren();
    }
}

template<int D>
void MRGrid<D>::saveGrid(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::loadGrid(const string &file) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRGrid<D>::yieldChildren(GridNodeVector &nodeTable, const NodeIndexSet &idxSet) {
    typename set<const NodeIndex<D> *>::iterator it;
    for (it = idxSet.begin(); it != idxSet.end(); it++) {
        GridNode<D> &node = this->rootBox->getNode(**it);
        int childDepth = node.getDepth() + 1;
        if (this->maxDepth != 0 and childDepth > this->maxDepth) {
            println(1, "+++ Maximum depth reached: " << childDepth);
            node.setIsEndNode();
        } else {
            node.createChildren();
            for (int i = 0; i < node.getNChildren(); i++) {
                GridNode<D> *child = node.getChild(i);
                nodeTable.push_back(child);
            }
        }
    }
}

template class MRGrid<1>;
template class MRGrid<2>;
template class MRGrid<3>;
