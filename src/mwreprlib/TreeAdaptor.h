#ifndef TREEADAPTOR_H
#define TREEADAPTOR_H

#include "mwrepr_declarations.h"

template<int D>
class TreeAdaptor {
public:
    TreeAdaptor() { }

    MRNodeVector* splitNodeVector(MRNodeVector &nodeVec,
                                  MRNodeVector *no_split = 0) const {
        MRNodeVector *split = new MRNodeVector;
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = static_cast<MWNode<D> &>(*nodeVec[n]);
            if (splitNode(node)) {
                split->push_back(&node);
            } else if (no_split != 0) {
                no_split->push_back(&node);
            }
        }
        return split;
    }

protected:
    virtual bool splitNode(MWNode<D> &node) const = 0;
};

#endif // TREEADAPTOR_H
