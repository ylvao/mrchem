#ifndef GRIDADAPTOR_H
#define GRIDADAPTOR_H

#include "mwrepr_declarations.h"

template<int D>
class GridAdaptor {
public:
    GridAdaptor(double pr = -1.0, bool abs = false);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    void adaptGrid(MRGrid<D> &g, FunctionTree<D> &t);

protected:
    bool absPrec;
    double prec;
    MRGrid<D> *grid;
    FunctionTree<D> *tree;

    void splitNodes(NodeIndexSet &splitSet, NodeIndexSet &gridSet);
    bool splitCheck(MWNode<D> &node);
    double getWaveletThreshold(int factor, int scale);
};

#endif // GRIDADAPTOR_H
