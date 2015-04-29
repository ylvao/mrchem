#ifndef HILBERTITERATOR_H
#define HILBERTITERATOR_H

#include "MRTree.h"
#include "MRNode.h"
#include "TreeIterator.h"
#include "HilbertPath.h"

template<int D>
class HilbertIterator: public TreeIterator<D> {
public:
    HilbertIterator(MRTree<D> *tree, int dir = TopDown):
        TreeIterator<D>(tree, dir) {
        this->init(tree);
    }
    virtual ~HilbertIterator() {}

protected:
    int getChildIndex(int i) const {
        const MRNode<D> &node = *this->state->node;
        const HilbertPath<D> &h = node.getHilbertPath();
        return h.getZIndex(i);
    }
};

#endif // HILBERTITERATOR_H
