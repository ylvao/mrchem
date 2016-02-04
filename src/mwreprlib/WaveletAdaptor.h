#ifndef WAVELETADAPTOR_H
#define WAVELETADAPTOR_H

#include "mwrepr_declarations.h"
#include "TreeAdaptor.h"

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr = -1.0, bool abs = false, int scale = 20)
            : prec(pr),
              absPrec(abs),
              maxScale(scale) { }

    virtual ~WaveletAdaptor() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setMaxScale(int scale) { this->maxScale = scale; }

protected:
    double prec;
    bool absPrec;
    int maxScale;

    virtual bool splitNode(MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // WAVELETADAPTOR_H
