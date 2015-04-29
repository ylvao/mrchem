#ifndef LEBESGUEITERATOR_H
#define LEBESGUEITERATOR_H

#include "MRTree.h"
#include "MRNode.h"
#include "TreeIterator.h"

template<int D>
class LebesgueIterator: public TreeIterator<D> {
public:
    LebesgueIterator(MRTree<D> *tree, int dir = TopDown):
        TreeIterator<D>(tree, dir) {
        this->init(tree);
    }
    virtual ~LebesgueIterator() {}
protected:
    int getChildIndex(int i) const { return i; }
};

#endif // LEBESGUEITERATOR_H
