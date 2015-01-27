#ifndef GRIDADAPTOR_H
#define GRIDADAPTOR_H

#include "mwrepr_declarations.h"

template<int D>
class GridAdaptor {
public:
    GridAdaptor(double pr = -1.0, bool abs = false);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    void adaptGrid(MRGrid<D> &grid, FunctionTree<D> &tree);
    MRNodeVector& splitNodeTable(MRNodeVector &nodeTable);

protected:
    bool absPrec;
    double prec;
};

#endif // GRIDADAPTOR_H
