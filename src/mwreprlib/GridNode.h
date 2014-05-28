/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#ifndef GRIDNODE_H_
#define GRIDNODE_H_

#include <Eigen/Core>

#include "macros.h"
#include "parallel.h"
#include "mwrepr_declarations.h"

#include "MRGrid.h"
#include "GridNodeBox.h"

#ifdef OPENMP
#define SET_NODE_LOCK() omp_set_lock(&this->node_lock)
#define UNSET_NODE_LOCK() omp_unset_lock(&this->node_lock)
#define TEST_NODE_LOCK() omp_test_lock(&this->node_lock)
#else
#define SET_NODE_LOCK()
#define UNSET_NODE_LOCK()
#define TEST_NODE_LOCK() false
#endif

template<int D>
class GridNode {
public:
    GridNode(MRGrid<D> *_grid, int n, const int *l);
    GridNode(GridNode<D> *_parent, const int *l);
    ~GridNode();

    void createChildren();
    void deleteChildren();

    void clearGridPointer() { this->grid = 0; }
    MRGrid<D> &getGrid() { return *this->grid; }

    GridNode<D> &getNode(const NodeIndex<D> &idx);

    int getTDim() const { return this->grid->getTDim(); }
    int getKp1() const { return this->grid->getKp1(); }
    int getKp1_d() const { return this->grid->getKp1_d(); }
    int getDepth() const { return this->nodeIndex.getScale()-this->grid->getRootScale(); }
    int getScale() const { return this->nodeIndex.getScale(); }
    int getOrder() const { return this->grid->getOrder(); }
    int getNChildren() const { if (isBranchNode()) return getTDim(); return 0; }
    void getCenter(double *r) const;
    void getLowerBounds(double *r) const;
    void getUpperBounds(double *r) const;
    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }
    const int *getTranslation() const { return this->nodeIndex.getTranslation(); }

    bool checkCoordOnNode(const double *r) const; 

    inline bool isRoot() const;
    inline bool hasCoefs() const;
    inline bool isEndNode() const;
    inline bool isGenNode() const;
    inline bool isLeafNode() const;
    inline bool isAllocated() const;
    inline bool isBranchNode() const;
    inline bool hasChild(int i) const;

    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }
    void setHasCoefs(bool flag = true) {
	if (flag) {
	    SET_BITS(status, FlagHasCoefs | FlagAllocated);
	} else {
	    CLEAR_BITS(status, FlagHasCoefs);
	}
    }

    const Eigen::MatrixXd &getQuadPoints() const { return this->roots; }
    const Eigen::MatrixXd &getQuadWeights() const { return this->weights; }

    void getExpandedPoints(Eigen::MatrixXd &points) const;
    void getExpandedWeights(Eigen::VectorXd &weights) const;

    GridNode<D> *getChild(int i) {
	assert(i >= 0 and i < getTDim());
	assert(this->children != 0);
	return this->children[i];
    }

    friend class MRGrid<D>;

protected:
    NodeIndex<D> nodeIndex;

    Eigen::MatrixXd roots;
    Eigen::MatrixXd weights;

    MRGrid<D> *grid;
    GridNode<D> *parent; ///< Parent node
    GridNode<D> **children; ///< 2^D children

    inline bool checkStatus(unsigned char mask) const;

    void createChild(int i);
    void allocKindergarten();

    GridNode<D> *retrieveNode(const NodeIndex<D> &idx);
    bool isAncestor(const NodeIndex<D> &idx) const;

    void calcChildTranslation(int cIdx, int *l) const;
    int getChildIndex(const NodeIndex<D> &nIdx) const;
    NodeIndex<D> getChildIndex(int cIdx) const;

    void calcQuadPoints();
    void calcQuadWeights();
//  void recurCompQuadRoots(int d, int *l, const Eigen::VectorXd &primRoots, Eigen::MatrixXd &expRoots);

    static const unsigned char FlagBranchNode =  B8(00000001);
    static const unsigned char FlagGenNode =     B8(00000010);
    static const unsigned char FlagHasCoefs =    B8(00000100);
    static const unsigned char FlagAllocated =   B8(00001000);
    static const unsigned char FlagEndNode =     B8(00010000);
    static const unsigned char FlagRootNode =    B8(00100000);
#ifdef OPENMP
    omp_lock_t node_lock;
#endif
private:
    unsigned char status;
};

template<int D>
bool GridNode<D>::isRoot() const {
    if (this->status & FlagRootNode) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::hasCoefs() const {
    if (this->status & FlagHasCoefs) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::isEndNode() const {
    if (this->status & FlagEndNode) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::isGenNode() const {
    if (this->status & FlagGenNode) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::isLeafNode() const {
    if (this->status & FlagBranchNode) {
	return false;
    }
    return true;
}

template<int D>
bool GridNode<D>::isAllocated() const {
    if (this->status & FlagAllocated) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::isBranchNode() const {
    if (this->status & FlagBranchNode) {
	return true;
    }
    return false;
}

template<int D>
bool GridNode<D>::hasChild(int i) const {
    assert(i >= 0 and i < getTDim());
    assert(this->children != 0);
    if (this->children[i] == 0) {
	return false;
    }
    return true;
}

template<int D>
bool GridNode<D>::checkStatus(unsigned char mask) const {
    if (mask == (this->status & mask)) {
    	return true;
    }
    return false;
}

#endif /* GRIDNODE_H_ */

