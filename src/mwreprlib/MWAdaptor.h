#ifndef MWADAPTOR_H
#define MWADAPTOR_H

#include "mwrepr_declarations.h"

template<int D>
class MWAdaptor {
public:
    MWAdaptor(double pr = -1.0, bool abs = false, int scale = 20);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    MRNodeVector* splitNodeVector(MRNodeVector &nodeVec,
                         MRNodeVector *noSplit = 0) const;

protected:
    int maxScale;
    bool absPrec;
    double prec;

    bool splitNode(MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // MWADAPTOR_H
