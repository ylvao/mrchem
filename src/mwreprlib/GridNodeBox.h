/**
 *
 *
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#ifndef GRIDNODEBOX_H_
#define GRIDNODEBOX_H_

#include "mwrepr_declarations.h"
#include "NodeIndex.h"

template<int D>
class GridNodeBox {
public:
    GridNodeBox(int scale = 0, const int *nbox = 0, const double *origo = 0);
    GridNodeBox(const GridNodeBox<D> &box);
    virtual ~GridNodeBox();

    void setNode(int idx, GridNode<D> *node);
    void removeNode(int idx);

    NodeIndex<D> getNodeIndex(const double *r);
    NodeIndex<D> getNodeIndex(int i);

    int getBoxIndex(const double *r);
    int getBoxIndex(const NodeIndex<D> &idx);
    
    GridNode<D> &getNode(const NodeIndex<D> &idx) const;
    GridNode<D> &getNode(const double *r) const;
    GridNode<D> &getNode(int i = 0) const;

    int getNBoxes(int d = -1) const;
    int getNOccupied() const { return this->nOccupied; }
    int getRootScale() const { return this->cornerIndex.getScale(); }
    double getUnitLength() const { return this->unitLength; }
    const double *getOrigin() const { return this->origin; }
    const double *getBoxLength() const { return this->boxLength; }
    const double *getLowerBounds() const { return this->lowerBounds; }
    const double *getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }
    GridNode<D> **getNodes() { return this->nodes; }

protected:
    int nBoxes[D+1];		///< Number of boxes in each dimension, last entry total
    int nOccupied;		///< Number of non-zero pointers in box
    double unitLength;		///< 1/2^initialScale
    double origin[D];		///< Relates box origin to real origin
    double boxLength[D];	///< Total length (unitLength times nBoxes)
    double lowerBounds[D];	///< Box lower bound (not real)
    double upperBounds[D];	///< Box upperBound (not real)
    GridNode<D> **nodes;	///< Container of root nodes
    NodeIndex<D> cornerIndex;	///< Index defining the lower corner of the box

    void allocNodePointers();
    void deleteNodes();
};

#endif /* GRIDNODEBOX_H_ */
