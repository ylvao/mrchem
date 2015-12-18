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
    BoundingBox(const NodeIndex<D> &idx, const int *nb = 0);
    virtual ~BoundingBox() { }

    BoundingBox(const BoundingBox<D> &box);
    BoundingBox<D> &operator=(const BoundingBox<D> &box);

    inline bool operator==(const BoundingBox<D> &box) const;
    inline bool operator!=(const BoundingBox<D> &box) const;

    NodeIndex<D> getNodeIndex(const double *r) const;
    NodeIndex<D> getNodeIndex(int bIdx) const;

    int getBoxIndex(const double *r) const;
    int getBoxIndex(const NodeIndex<D> &nIdx) const;

    int getNBoxes() const { return this->nBoxes[D]; }
    int getNBoxes(int d) const { return this->nBoxes[d]; }
    int getRootScale() const { return this->cornerIndex.getScale(); }
    double getUnitLength() const { return this->unitLength; }
    double getBoxLength(int d) const { return this->boxLengths[d]; }
    double getLowerBound(int d) const { return this->lowerBounds[d]; }
    double getUpperBound(int d) const { return this->upperBounds[d]; }
    const double *getBoxLengths() const { return this->boxLengths; }
    const double *getLowerBounds() const { return this->lowerBounds; }
    const double *getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const BoundingBox<T> &box);

protected:
    // Fundamental parameters
    int nBoxes[D+1];		        ///< Number of boxes in each dim, last entry total
    NodeIndex<D> cornerIndex;	///< Index defining the lower corner of the box

    // Derived parameters
    double unitLength;		    ///< 1/2^initialScale
    double boxLengths[D];	        ///< Total length (unitLength times nBoxes)
    double lowerBounds[D];	    ///< Box lower bound (not real)
    double upperBounds[D];	    ///< Box upper bound (not real)

    void setNBoxes(const int *nb);
    void setDerivedParameters();

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & nBoxes;
        ar & cornerIndex;
        ar & unitLength;
        ar & boxLengths;
        ar & lowerBounds;
        ar & upperBounds;
    }
};

template<int D>
bool BoundingBox<D>::operator==(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return false;
    for (int d = 0; d < 3; d++) {
        if (getNBoxes(d) != box.getNBoxes(d)) return false;
    }
    return true;
}

template<int D>
bool BoundingBox<D>::operator!=(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return true;
    for (int d = 0; d < 3; d++) {
        if (getNBoxes(d) != box.getNBoxes(d)) return true;
    }
    return false;
}

template<int T>
std::ostream& operator<<(std::ostream &o, const BoundingBox<T> &box) {
    GET_PRINT_PRECISION(int pprec);
    o << std::fixed << std::setprecision(5);
    o << "*BoundingBox: " << std::endl;
    o << "  unit length     = " << box.getUnitLength() << std::endl;
    o << "  total boxes     = " << box.getNBoxes() << std::endl;
    o << "  boxes           = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(11) << box.getNBoxes(i) << " ";
    }
    o << " ]" << std::endl;
    o << "  lower bounds    = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(11) << box.getLowerBound(i) << " ";
    }
    o << " ]" << std::endl;
    o << "  upper bounds    = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(11) << box.getUpperBound(i) << " ";
    }
    o << " ]" << std::endl;
    o << "  total length    = [ ";
    for (int i = 0; i < T; i++) {
        o << std::setw(11) << box.getBoxLength(i) << " ";
    }
    o << " ]" << std::endl;
    o << std::scientific << std::setprecision(pprec);
    return o;
}

#endif /* BOUNDINGBOX_H_ */
