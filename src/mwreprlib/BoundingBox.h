/**
 *
 *
 *          CTCC, University of Tromsø
 *
 */

#ifndef BOUNDINGBOX_H_
#define BOUNDINGBOX_H_

#include "mwrepr_declarations.h"
#include "NodeIndex.h"

template<int D>
class BoundingBox {
public:
    BoundingBox(const NodeIndex<D> &idx, const int *nb = 0, const double *o = 0);
    virtual ~BoundingBox() { }

    BoundingBox(const BoundingBox<D> &box);
    BoundingBox<D> &operator=(const BoundingBox<D> &box);

    bool operator==(const BoundingBox<D> &box) const { NOT_IMPLEMENTED_ABORT; }
    bool operator!=(const BoundingBox<D> &box) const { NOT_IMPLEMENTED_ABORT; }

    void setCornerIndex(const NodeIndex<D> &idx);
    void setOrigin(const double *o);
    void setNBoxes(const int *nb);

    NodeIndex<D> getNodeIndex(const double *r) const;
    NodeIndex<D> getNodeIndex(int bIdx) const;

    int getBoxIndex(const double *r) const;
    int getBoxIndex(const NodeIndex<D> &nIdx) const;

    inline int getNBoxes(int d = -1) const;
    int getRootScale() const { return this->cornerIndex.getScale(); }
    double getUnitLength() const { return this->unitLength; }
    const double *getOrigin() const { return this->origin; }
    const double *getBoxLength() const { return this->boxLength; }
    const double *getLowerBounds() const { return this->lowerBounds; }
    const double *getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const BoundingBox<T> &box);

protected:
    // Fundamental parameters
    int nBoxes[D+1];		///< Number of boxes in each dim, last entry total
    double origin[D];		///< Relates box origin to real origin
    NodeIndex<D> cornerIndex;	///< Index defining the lower corner of the box

    // Derived parameters
    double unitLength;		///< 1/2^initialScale
    double boxLength[D];	///< Total length (unitLength times nBoxes)
    double lowerBounds[D];	///< Box lower bound (not real)
    double upperBounds[D];	///< Box upper bound (not real)

    void setDerivedParameters();

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & nBoxes;
        ar & unitLength;
        ar & origin;
        ar & boxLength;
        ar & lowerBounds;
        ar & upperBounds;
        ar & cornerIndex;
    }
};


template<int D>
int BoundingBox<D>::getNBoxes(int d) const {
    if (d < 0) {
        return this->nBoxes[D];
    } else if (d < D) {
        return this->nBoxes[d];
    } else {
        MSG_ERROR("Invalid dimension argument");
        return -1;
    }
}

template<int T>
std::ostream& operator<<(std::ostream &o, const BoundingBox<T> &box) {
    GET_PRINT_PRECISION(int pprec);
    o << std::fixed;
    o << "*BoundingBox: " << std::endl;
    o << "  corner index    = " << box.cornerIndex << std::endl;
    o << "  unit box length = " << box.unitLength << std::endl;
    o << "  origin          = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(21) << box.origin[i] << " ";
    }
    o << "]" << std::endl;
    o << "  lower bounds    = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(21) << box.lowerBounds[i] << " ";
    }
    o << "]" << std::endl;
    o << "  upper bounds    = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(21) << box.upperBounds[i] << " ";
    }
    o << "]" << std::endl;
    o << "  boxes           = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(21) << box.nBoxes[i] << " ";
    }
    o << "]" << std::endl;
    o << "  nBoxes          = " << box.nBoxes[T] << std::endl;
    o << std::scientific << std::setprecision(pprec);
    return o;
}

#endif /* BOUNDINGBOX_H_ */
