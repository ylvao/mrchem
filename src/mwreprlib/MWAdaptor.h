#ifndef MWADAPTOR_H
#define MWADAPTOR_H

#include "mwrepr_declarations.h"

template<int D>
class MWAdaptor {
public:
    MWAdaptor(double pr = -1.0, bool abs = false, int scale = 20);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    void splitNodeVector(double norm, MRNodeVector &nodeVec,
                                      MRNodeVector &splitVec,
                                      MRNodeVector &noSplitVec);

protected:
    int maxScale;
    bool absPrec;
    double prec;

    bool splitCheck(double norm, MWNode<D> &node);
    double getWaveletThreshold(double norm, int scale);
};

#endif // MWADAPTOR_H
