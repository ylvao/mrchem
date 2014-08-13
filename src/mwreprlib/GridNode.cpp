/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include "MathUtils.h"
#include "GridNode.h"
#include "NodeIndex.h"
#include "QuadratureCache.h"

using namespace Eigen;
using namespace std;

template<int D>
GridNode<D>::GridNode(MRTree<D> *t, const NodeIndex<D> &idx)
    : MRNode<D>(t, idx) {
    calcQuadPoints();
    calcQuadWeights();
}

template<int D>
GridNode<D>::GridNode(GridNode<D> *p, const int *l) : MRNode<D>(p, l) {
    calcQuadPoints();
    calcQuadWeights();
}

template<int D>
GridNode<D>::~GridNode() {
    if (this->isGenNode()) {
        this->tree->decrementGenNodeCount();
        this->tree->decrementAllocGenNodeCount();
    }
}

template<int D>
void GridNode<D>::createChild(int i) {
    assert(this->children[i] == 0);
    int l[3];
    this->calcChildTranslation(i, l);
    GridNode<D> *child = new GridNode<D> (this, l);
    this->children[i] = child;
    child->setIsEndNode();
}

template<int D>
void GridNode<D>::genChild(int i) {
    assert(this->children[i] == 0);
    int l[3];
    this->calcChildTranslation(i, l);
    GridNode<D> *child = new GridNode<D> (this, l);
    this->children[i] = child;
    child->setIsGenNode();
    this->tree->incrementGenNodeCount();
    this->tree->incrementAllocGenNodeCount();
}

template<int D>
MRNode<D> *GridNode<D>::retrieveNode(int n, const double *r) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D> *GridNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (this->nodeIndex.getScale() == idx.getScale()) { // we're done
        return this;
    }
    SET_NODE_LOCK();
    if (this->isLeafNode()) {
        this->genChildren();
    }
    UNSET_NODE_LOCK();
    int cIdx = this->getChildIndex(idx);
    return this->children[cIdx]->retrieveNode(idx);
}

template<int D>
void GridNode<D>::calcQuadPoints() {
    int kp1 = this->getKp1();
    this->roots = MatrixXd::Zero(kp1,D);

    getQuadratureCache(qc);
    const VectorXd &pts = qc.getRoots(kp1);

    double sFac = pow(2.0, -this->getScale());
    const int *l = this->getTranslation();
    const double *o = this->tree->getOrigin();
    for (int d = 0; d < D; d++) {
        this->roots.col(d) = sFac*(pts.array() + double(l[d])) - o[d];
    }
}

template<int D>
void GridNode<D>::calcQuadWeights() {
    int kp1 = this->getKp1();
    double sFac = pow(2.0, -this->getScale());
    this->weights = MatrixXd::Zero(kp1,D);

    getQuadratureCache(qc);
    VectorXd wgts = sFac*qc.getWeights(kp1);

    for (int d = 0; d < D; d++) {
        this->weights.col(d) = wgts;
    }
}

template<int D>
void GridNode<D>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    NOT_IMPLEMENTED_ABORT;
}

template<>
void GridNode<1>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    int kp1_d = this->getKp1_d();
    expandedPoints = MatrixXd::Zero(kp1_d,1);

    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints.col(0) = primitivePoints;
}

template<>
void GridNode<2>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    NOT_IMPLEMENTED_ABORT;
    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints = MatrixXd::Zero(kp1_d,2);
    MathUtils::tensorExpandCoords_2D(kp1, primitivePoints, expandedPoints);
}

template<>
void GridNode<3>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints = MatrixXd::Zero(kp1_d,3);
    MathUtils::tensorExpandCoords_3D(kp1, primitivePoints, expandedPoints);
}

template<int D>
void GridNode<D>::getExpandedWeights(Eigen::VectorXd &expandedWeights) const {
    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    int inpos = kp1_d - kp1;
    const MatrixXd &primitiveWeights = getQuadWeights();

    expandedWeights = VectorXd::Zero(kp1_d);
    for (int i = 0; i < kp1; i++) {
        expandedWeights(inpos + i) = primitiveWeights(i,0);
    }
    MathUtils::tensorExpandCoefs(D, 0, kp1, kp1_d, primitiveWeights, expandedWeights);
}

template class GridNode<1>;
template class GridNode<2>;
template class GridNode<3>;
