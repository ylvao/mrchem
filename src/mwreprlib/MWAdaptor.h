#ifndef MWADAPTOR_H
#define MWADAPTOR_H

#include "mwrepr_declarations.h"
#include "TreeAdaptor.h"

template<int D>
class MWAdaptor : public TreeAdaptor<D> {
public:
    MWAdaptor(double pr = -1.0, bool abs = false, int scale = 20);

    void setMaxScale(int scale) { this->maxScale = scale; }
    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

protected:
    int maxScale;
    bool absPrec;
    double prec;

    bool splitNode(MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // MWADAPTOR_H
