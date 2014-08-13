/**
 *
 *
 *          CTCC, University of Tromsø
 *
 */

#include "MathUtils.h"
#include "constants.h"
#include "BoundingBox.h"

using namespace std;

template<int D>
BoundingBox<D>::BoundingBox(const NodeIndex<D> &idx, const int *nb,
                            const double *o) : cornerIndex(idx) {
    this->nBoxes[D] = 1;
    for (int d = 0; d < D; d++) {
        if (nb == 0) {
            this->nBoxes[d] = 1;
        } else {
            assert(nb[d] > 0);
            this->nBoxes[d] = nb[d];
            this->nBoxes[D] *= this->nBoxes[d];
        }
        if (o == 0) {
            this->origin[d] = 0.0;
        } else {
            this->origin[d] = o[d];
        }
    }
    setDerivedParameters();
}

template<int D>
BoundingBox<D>::BoundingBox(const BoundingBox<D> &box)
        : cornerIndex(box.cornerIndex) {
    this->nBoxes[D] = box.nBoxes[D];
    for (int d = 0; d < D; d++) {
        assert(box.nBoxes[d] > 0);
        this->nBoxes[d] = box.nBoxes[d];
        this->origin[d] = box.origin[d];
    }
    setDerivedParameters();
}

template<int D>
BoundingBox<D> &BoundingBox<D>::operator=(const BoundingBox<D> &box) {
    this->cornerIndex = box.cornerIndex;
    this->nBoxes[D] = box.nBoxes[D];
    for (int d = 0; d < D; d++) {
        assert(box.nBoxes[d] > 0);
        this->nBoxes[d] = box.nBoxes[d];
        this->origin[d] = box.origin[d];
    }
    setDerivedParameters();
}

template<int D>
void BoundingBox<D>::setCornerIndex(const NodeIndex<D> &idx) {
    this->cornerIndex = idx;
    setDerivedParameters();
}

template<int D>
void BoundingBox<D>::setOrigin(const int *o) {
    for (int d = 0; d < D; d++) {
        if (o == 0) {
            this->origin[d] = 0.0;
        } else {
            this->origin[d] = o[d];
        }
    }
    setDerivedParameters();
}

template<int D>
void BoundingBox<D>::setNBoxes(const int *nb) {
    this->nBoxes[D] = 1;
    for (int d = 0; d < D; d++) {
        if (nb == 0) {
            this->nBoxes[d] = 1;
        } else {
            assert(nb[d] > 0);
            this->nBoxes[d] = nb[d];
            this->nBoxes[D] *= this->nBoxes[d];
        }
    }
    setDerivedParameters();
}

template<int D>
void BoundingBox<D>::setDerivedParameters() {
    assert(this->nBoxes[D] > 0);
    int scale = this->cornerIndex.getScale();
    const int *l = this->cornerIndex.getTranslation();
    this->unitLength = pow(2.0, -scale);
    for (int d = 0; d < D; d++) {
        assert(box.nBoxes[d] > 0);
        this->boxLength[d] = this->unitLength * this->nBoxes[d];
        this->lowerBounds[d] = l[d] * this->unitLength;
        this->upperBounds[d] = this->lowerBounds[d] + this->boxLength[d];
    }
}

template<int D>
NodeIndex<D> BoundingBox<D>::getNodeIndex(const double *r) const {
    assert(r != 0);
    int idx[D];
    for (int d = 0; d < D; d++) {
        double x = r[d];
        assert(x >= this->lowerBounds[d]);
        assert(x < this->upperBounds[d]);
        double div = (x - this->lowerBounds[d]) / this->unitLength;
        double iint;
        double ffloat = modf(div,&iint);
        idx[d] = (int) iint;
    }

    const int *cl = this->cornerIndex.getTranslation();

    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = idx[d] + cl[d];
    }

    int n = getRootScale();
    NodeIndex<D> nIdx(n, l);
    return nIdx;
}

template<>
NodeIndex<1> BoundingBox<1>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= nBoxes[1]);
    int n = getRootScale();
    int cl = this->cornerIndex.getTranslation(0);
    int l = bIdx + cl;
    NodeIndex<1> nIdx(n, &l);
    return nIdx;
}

template<int D>
NodeIndex<D> BoundingBox<D>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= nBoxes[D]);
    int l[D];
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) {
            ncells *= this->nBoxes[i];
        }
        double div = bIdx / ncells;
        double iint;
        modf(div, &iint);
        l[d] = (int) iint;
        bIdx -= ncells * l[d];
    }

    int n = getRootScale();
    const int *cl = this->cornerIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        l[d] += cl[d];
    }
    NodeIndex<D> nIdx(n, l);
    return nIdx;
}

template<>
int BoundingBox<1>::getBoxIndex(const double *r) const {
    assert(r != 0);
    double x = r[0];
    assert(x >= this->lowerBounds[0]);
    assert(x < this->upperBounds[0]);
    double div = (x - this->lowerBounds[0]) / this->unitLength;
    double iint;
    double ffloat = modf(div,&iint);
    return (int) iint;
}

template<int D>
int BoundingBox<D>::getBoxIndex(const double *r) const {
    assert(r != 0);
    int idx[D];
    for (int d = 0; d < D; d++) {
        double x = r[d];
        assert(x >= this->lowerBounds[d]);
        assert(x < this->upperBounds[d]);
        double div = (x - this->lowerBounds[d]) / this->unitLength;
        double iint;
        double ffloat = modf(div,&iint);
        idx[d] = (int) iint;
    }

    int bIdx = 0;
    for (int i = D - 1; i >= 0; i--) {
        int ncells = 1;
        for (int j = 0; j < i; j++) {
            ncells *= nBoxes[j];
        }
        bIdx += ncells * idx[i];
    }
    return bIdx;
}

template<>
int BoundingBox<1>::getBoxIndex(const NodeIndex<1> &nIdx) const {
    int n = nIdx.getScale();
    int l = nIdx.getTranslation(0);
    int cn = this->cornerIndex.getScale();
    int cl = this->cornerIndex.getTranslation(0);
    int relScale = n - cn;
    assert(relScale >= 0);
    return (l >> relScale) - cl;
}

template<int D>
int BoundingBox<D>::getBoxIndex(const NodeIndex<D> &nIdx) const {
    int n = nIdx.getScale();
    int cn = this->cornerIndex.getScale();
    const int *l = nIdx.getTranslation();
    const int *cl = this->cornerIndex.getTranslation();
    int relScale = n - cn;
    assert(relScale >= 0);

    int bIdx = 0;
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) {
            ncells *= this->nBoxes[i];
        }
        int reqTransl = (l[d] >> relScale) - cl[d];
        bIdx += ncells * reqTransl;
    }
    return bIdx;
}

template class BoundingBox<1>;
template class BoundingBox<2>;
template class BoundingBox<3>;
