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

protected:
    Eigen::MatrixXd roots;
    Eigen::MatrixXd weights;

    MRNode<D> *retrieveNode(int n, const double *r);
    MRNode<D> *retrieveNode(const NodeIndex<D> &idx);

    void createChild(int cIdx);
    void genChild(int cIdx);

    void calcQuadPoints();
    void calcQuadWeights();
};

#endif /* GRIDNODE_H_ */

