#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "mwrepr_declarations.h"

template<int D>
class MWProjector {
public:
    MWProjector();
    MWProjector(MWAdaptor<D> &a);
    virtual ~MWProjector();

protected:
    MWTree<D> *outTree;
    MWAdaptor<D> *adaptor;

    void buildTree();
    void calcNodeVector(MRNodeVector &nodeVec);
    virtual void calcWaveletCoefs(MWNode<D> &node) = 0;

    MRNodeVector* clearForeignNodes(MRNodeVector *oldVec) const;
    NodeIndexSet* getNodeIndexSet(const MRNodeVector &nodeVec) const;
};

#endif // MWPROJECTOR_H
