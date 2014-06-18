/**
 *
 *
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#ifndef NODEBOX_H_
#define NODEBOX_H_

#include "BoundingBox.h"

template<int D>
class NodeBox : public BoundingBox<D> {
public:
    NodeBox(int scale = 0, const int *nbox = 0, const double *origo = 0);
    NodeBox(const NodeBox<D> &box);
    virtual ~NodeBox();

    void initBox(const int *nbox, const double *origo);

    void setNode(int idx, MRNode<D> **node);
    void removeNode(int idx);

    MRNode<D> &getNode(const NodeIndex<D> &idx);
    MRNode<D> &getNode(const double *r);
    MRNode<D> &getNode(int i = 0);

    int getNOccupied() const { return this->nOccupied; }
    MRNode<D> **getNodes() { return this->nodes; }

protected:
    int nOccupied;		///< Number of non-zero pointers in box
    MRNode<D> **nodes;		///< Container of nodes

    void allocNodePointers();
    void deleteNodes();
};

#endif /* NODEBOX_H_ */
