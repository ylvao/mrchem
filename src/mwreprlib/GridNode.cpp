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
GridNode<D>::GridNode(MRGrid<D> *_grid, int n, const int *l) {
    if (_grid == 0) {
        MSG_FATAL("Cannot initialize node without tree!");
    }
    this->nodeIndex.setScale(n);
    this->nodeIndex.setTranslation(l);

    this->grid = _grid;
    this->grid->incrementNodeCount(n);

    this->parent = 0;
    this->status = 0;
    this->children = 0;

    calcQuadPoints();
    calcQuadWeights();

    setIsLeafNode();
    setIsEndNode();
    setIsRootNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
GridNode<D>::GridNode(GridNode<D> *_parent, const int *l) {
    this->parent = _parent;
    this->status = 0;
    this->children = 0;

    if (this->parent == 0) {
	NOT_IMPLEMENTED_ABORT
    } else {
	int n = this->parent->getScale() + 1;
        this->nodeIndex.setScale(n);
	this->nodeIndex.setTranslation(l);
        this->grid = this->parent->grid;
	this->grid->incrementNodeCount(n);
    }

    calcQuadPoints();
    calcQuadWeights();

    setIsLeafNode();
    setIsEndNode();
#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
GridNode<D>::~GridNode() {
    SET_NODE_LOCK();
    if (this->isBranchNode()) {
        assert(this->children != 0);
        deleteChildren();
    }
    this->grid->decrementNodeCount(getScale());
    UNSET_NODE_LOCK();
#ifdef OPENMP
    omp_destroy_lock(&node_lock);
#endif
}

template<int D>
void GridNode<D>::deleteChildren() {
    assert(this->children != 0);
    for (int n = 0; n < getTDim(); n++) {
        if (this->children[n] != 0) {
            delete this->children[n];
            this->children[n] = 0;
        }
    }
    delete [] this->children;
    this->children = 0;
    this->setIsLeafNode();
}


template<int D>
void GridNode<D>::createChildren() {
    if (this->children == 0) {
        this->allocKindergarten();
    }
    for (int n = 0; n < getTDim(); n++) {
        createChild(n);
    }
    this->setIsBranchNode();
    this->clearIsEndNode();
}

template<int D>
void GridNode<D>::createChild(int i) {
    assert(this->children[i] == 0);
    int l[3];
    calcChildTranslation(i, l);
    GridNode<D> *child = new GridNode<D> (this, l);
    this->children[i] = child;
    child->setIsEndNode();
}

template<int D>
void GridNode<D>::allocKindergarten() {
    if (this->children == 0) {
	int nChildren = getTDim();
	this->children = new GridNode<D> *[nChildren];
	for (int n = 0; n < nChildren; n++) {
	    this->children[n] = 0;
	}
    }
}

template<int D>
void GridNode<D>::calcQuadPoints() {
    int kp1 = getKp1();
    this->roots = MatrixXd::Zero(kp1,D);

    getQuadratureCache(qc);
    const VectorXd &pts = qc.getRoots(kp1);

    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    const double *o = this->grid->getOrigin();
    for (int d = 0; d < D; d++) {
	this->roots.col(d) = sFac*(pts.array() + double(l[d])) - o[d];
    }
}

template<int D>
void GridNode<D>::calcQuadWeights() {
    int kp1 = getKp1();
    double sFac = pow(2.0, -getScale());
    this->weights = MatrixXd::Zero(kp1,D);

    getQuadratureCache(qc);
    VectorXd wgts = sFac*qc.getWeights(kp1);

    for (int d = 0; d < D; d++) {
        this->weights.col(d) = wgts;
    }
}

template<int D>
void GridNode<D>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    NOT_IMPLEMENTED_ABORT
}

template<>
void GridNode<1>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    int kp1_d = getKp1_d();
    expandedPoints = MatrixXd::Zero(kp1_d,1);

    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints.col(0) = primitivePoints;
}

template<>
void GridNode<2>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    NOT_IMPLEMENTED_ABORT
    int kp1 = getKp1();
    int kp1_d = getKp1_d();
    int inpos = kp1_d - kp1;
    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints = MatrixXd::Zero(kp1_d,2);
    MathUtils::tensorExpandCoords_2D(kp1, primitivePoints, expandedPoints);
}

template<>
void GridNode<3>::getExpandedPoints(Eigen::MatrixXd &expandedPoints) const {
    int kp1 = getKp1();
    int kp1_d = getKp1_d();
    int inpos = kp1_d - kp1;
    const MatrixXd &primitivePoints = getQuadPoints();
    expandedPoints = MatrixXd::Zero(kp1_d,3);
    MathUtils::tensorExpandCoords_3D(kp1, primitivePoints, expandedPoints);
}

template<int D>
void GridNode<D>::getExpandedWeights(Eigen::VectorXd &expandedWeights) const {
    int kp1 = getKp1();
    int kp1_d = getKp1_d();
    int inpos = kp1_d - kp1;
    const MatrixXd &primitiveWeights = getQuadWeights();

    expandedWeights = VectorXd::Zero(kp1_d);
    expandedWeights.segment(inpos, kp1) = primitiveWeights.col(0);
    MathUtils::tensorExpandCoefs(D, 0, kp1, kp1_d, primitiveWeights, expandedWeights);
}

template<int D>
GridNode<D> &GridNode<D>::getNode(const NodeIndex<D> &idx) {
    if (not isAncestor(idx)) {
	THROW_ERROR("Cannot get node, node out of bounds:\nthis->idx: " <<
                        this->getNodeIndex() << "\narg->idx: " << idx);
    }
    int scale = idx.scale();
    if (scale > this->grid->getMaxScale()) {
	THROW_ERROR("Requested (" << scale << ") node beyond max scale ("
                        << this->grid->getMaxScale() << ")!")
    } 
    return *(retrieveNode(idx));
}

template<int D>
GridNode<D> *GridNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (this->nodeIndex.getScale() == idx.getScale()) { // we're done
        return this;
    }
    if (this->isLeafNode()) {
	THROW_ERROR("Requested node does not exist");
    }
    int childIdx = getChildIndex(idx);
    return children[childIdx]->retrieveNode(idx);
}

template<int D>
bool GridNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int inpScale = idx.scale();
    int relScale = inpScale - this->nodeIndex.scale();
    if (relScale < 0) {
        return false;
    }
    for (int d = 0; d < D; d++) {
        NodeIndex<D> &l = const_cast<NodeIndex<D> &> (idx);
        int reqTransl = l[d] >> relScale;
	if (this->nodeIndex[d] != reqTransl) {
	    return false;
	}
    }
    return true;
}

template<int D>
int GridNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    int cIdx = 0;
    int delta_scale = nIdx.getScale() - this->nodeIndex.getScale() - 1;
    for (int d = 0; d < D; d++) {
        int bit = (nIdx.getTranslation()[d] >> (delta_scale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    return cIdx;
}

template<int D>
NodeIndex<D> GridNode<D>::getChildIndex(int cIdx) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GridNode<D>::calcChildTranslation(int cIdx, int *l) const {
    for (int d = 0; d < D; d++) {
        l[d] = (2 * this->nodeIndex[d]) + ((cIdx >> d) & 1);
    }
}

template<int D>
void GridNode<D>::getCenter(double *r) const {
    if (r == 0) {
	THROW_ERROR("Invalid argument");
    }
    const double *origin = this->grid->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
	r[d] = sFac*(1.0*l[d] + 0.5) - origin[d];
    }
}

template<int D>
void GridNode<D>::getLowerBounds(double *r) const {
    if (r == 0) {
	THROW_ERROR("Invalid argument");
    }
    const double *origin = this->grid->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
	r[d] = sFac*(1.0*l[d]) - origin[d];
    }
}

template<int D>
void GridNode<D>::getUpperBounds(double *r) const {
    if (r == 0) {
	THROW_ERROR("Invalid argument");
    }
    const double *origin = this->grid->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
	r[d] = sFac*(1.0*l[d] + 1.0) - origin[d];
    }
}

template<int D>
bool GridNode<D>::checkCoordOnNode(const double *r) const {
    int n = getScale();
    double nLength = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
	const double *origin = this->grid->getOrigin();
	if (r[d] < (l[d]*nLength - origin[d]) or r[d] > ((l[d] + 1)*nLength) - origin[d]) {
            return false;
        }
    }
    return true;
}

template class GridNode<1>;
template class GridNode<2>;
template class GridNode<3>;
