#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

template<int D> class MRGrid;
template<int D> class GridAdaptor;

template<int D>
class MWProjector {
public:
    MWProjector(MRGrid<D> *startGrid, GridAdaptor<D> *adap);

protected:
    MRGrid<D> *grid;
    GridAdaptor<D> *adaptor;
};

#endif // MWPROJECTOR_H
