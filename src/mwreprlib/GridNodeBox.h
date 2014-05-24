/**
 *
 *
 *  \date Aug 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif This class is a container of Nodes defined in a region of space.
 * The NodeBox has tow operational modes:
 *   1) When the tree pointer is set, the box acts as the owner of the nodes
 *       stored in it (i.e. root nodes), and when the box is destroyed the
 *       nodes are destructed as well.
 *   2) When the tree pointer is null, the NodeBox acts just like a container
 *
 * The implementation tries to be computationally efficient, since these boxes
 * are created, used, and destroyed extensively during operator application.
 */

#ifndef GRIDNODEBOX_H_
#define GRIDNODEBOX_H_

#include "mwrepr_declarations.h"

template<int D>
class GridNodeBox {
public:
    GridNodeBox(int scale = 0, const int *nbox = 0, const double *origo = 0);
    virtual ~GridNodeBox();

    GridNode<D> &getNode(const NodeIndex<D> &idx) const;
    GridNode<D> &getNode(const double *r) const;
    GridNode<D> &getNode(int i = 0) const;

    void setNode(int idx, GridNode<D> *node);
    void removeNode(int idx);

    int getNBoxes(int d = -1) const;
    int getNOccupied() const { return this->nOccupied; }
    double getUnitLength() const { return this->unitLength; }
    const double *getOrigin() const { return this->origin; }
    const double *getBoxLength() const { return this->boxLength; }
    const double *getLowerBounds() const { return this->lowerBounds; }
    const double *getUpperBounds() const { return this->upperBounds; }
    NodeIndex<D> getNodeIndex() const { return this->cornerIndex; }
    
protected:
    int nBoxes[D+1];		///< Number of boxes in each dimension, last entry total
    int nOccupied;		///< Number of non-zero pointers in box
    double unitLength;		///< 1/2^initialScale
    double origin[D];		///< Relates box origin to real origin
    double boxLength[D];	///< Total length (unitLength times nBoxes)
    double lowerBounds[D];	///< Box lower bound (not real)
    double upperBounds[D];	///< Box upperBound (not real)
    GridNode<D> **nodes;	///< Container of root nodes
    NodeIndex<D> cornerIndex;

    void allocNodePointers(int n = 0);
    void deleteNodes();
};

#endif /* GRIDNODEBOX_H_ */
