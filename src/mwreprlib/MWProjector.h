#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

template<int D> class MRGrid;
template<int D> class GridAdaptor;
template<int D> class FunctionTree;

template<int D>
class MWProjector {
public:
    MWProjector();
    MWProjector(GridAdaptor<D> &a);
    virtual ~MWProjector();

protected:
    GridAdaptor<D> *adaptor;

    void buildTree(FunctionTree<D> &tree);
};

#endif // MWPROJECTOR_H
