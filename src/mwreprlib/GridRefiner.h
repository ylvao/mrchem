#ifndef GRIDADAPTOR_H
#define GRIDADAPTOR_H

#include "GridGenerator.h"

template<int D> class FunctionTree;

template<int D>
class GridAdaptor : public GridGenerator<D> {
public:
    GridAdaptor(double pr, bool abs = false);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    void adaptGrid(MRGrid<D> &outGrid, FunctionTree<D> &tree);

protected:
    bool absPrec;
    double prec;

    bool splitCheck(const MRNode<D> *node);
};

#endif // GRIDADAPTOR_H
