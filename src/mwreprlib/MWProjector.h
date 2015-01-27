#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "mwrepr_declarations.h"

template<int D>
class MWProjector {
public:
    MWProjector();
    MWProjector(GridAdaptor<D> &a);
    virtual ~MWProjector();

protected:
    FunctionTree<D> *tree;
    GridAdaptor<D> *adaptor;

    void buildTree();
    void calcNodeTable(MRNodeVector &nodeTable);
    virtual void calcWaveletCoefs(MWNode<D> &node) = 0;
};

#endif // MWPROJECTOR_H
