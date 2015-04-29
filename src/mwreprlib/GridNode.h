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

#include "MRNode.h"
#include "MRGrid.h"

template<int D>
class GridNode : public MRNode<D> {
public:
    GridNode(MRTree<D> &t, const NodeIndex<D> &idx);
    GridNode(GridNode<D> *p, int cIdx);
    virtual ~GridNode();

    const Eigen::MatrixXd &getQuadPoints() const { return this->roots; }
    const Eigen::MatrixXd &getQuadWeights() const { return this->weights; }

    void getExpandedPoints(Eigen::MatrixXd &points) const;
    void getExpandedWeights(Eigen::VectorXd &weights) const;

    inline MRGrid<D> &getGridTree();
    inline GridNode<D> &getGridParent();
    inline GridNode<D> &getGridChild(int cIdx);

    inline const MRGrid<D> &getGridTree() const;
    inline const GridNode<D> &getGridParent() const;
    inline const GridNode<D> &getGridChild(int cIdx) const;

protected:
    Eigen::MatrixXd roots;
    Eigen::MatrixXd weights;

    MRNode<D> *retrieveNode(int n, const double *r);
    MRNode<D> *retrieveNode(const NodeIndex<D> &idx);

    void createChild(int cIdx);
    void genChild(int cIdx);

    void calcQuadPoints();
    void calcQuadWeights();

    mpi::request isendCoefs(int who, int tag, int comp = -1);
    mpi::request ireceiveCoefs(int who, int tag, int comp = -1);
};

template<int D>
const MRGrid<D> &GridNode<D>::getGridTree() const {
    assert(this->tree != 0);
    return static_cast<const MRGrid<D> &>(*this->tree);
}

template<int D>
MRGrid<D> &GridNode<D>::getGridTree() {
    assert(this->tree != 0);
    return static_cast<MRGrid<D> &>(*this->tree);
}

template<int D>
const GridNode<D> &GridNode<D>::getGridChild(int cIdx) const {
    assert(this->children != 0);
    assert(this->children[cIdx] != 0);
    return static_cast<const GridNode<D> &>(*this->children[cIdx]);
}

template<int D>
GridNode<D> &GridNode<D>::getGridChild(int cIdx) {
    assert(this->children != 0);
    assert(this->children[cIdx] != 0);
    return static_cast<GridNode<D> &>(*this->children[cIdx]);
}

template<int D>
const GridNode<D> &GridNode<D>::getGridParent() const {
    assert(this->parent != 0);
    return static_cast<const GridNode<D> &>(*this->parent);
}

template<int D>
GridNode<D> &GridNode<D>::getGridParent() {
    assert(this->parent != 0);
    return static_cast<GridNode<D> &>(*this->parent);
}

#endif /* GRIDNODE_H_ */

