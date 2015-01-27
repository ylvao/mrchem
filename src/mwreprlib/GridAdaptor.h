#ifndef GRIDADAPTOR_H
#define GRIDADAPTOR_H

template<int D> class MRGrid;
template<int D> class FunctionTree;

template<int D>
class GridAdaptor {
public:
    GridAdaptor(double pr = -1.0, bool abs = false);

    void setAbsPrec(bool abs) { this->absPrec = abs; }
    void setPrecision(double pr) { this->prec = pr; }

    void adaptGrid(MRGrid<D> &grid, FunctionTree<D> &tree);

protected:
    bool absPrec;
    double prec;
};

#endif // GRIDADAPTOR_H
