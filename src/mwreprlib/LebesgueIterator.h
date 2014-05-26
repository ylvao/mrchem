#ifndef LEBESGUEITERATOR_H
#define LEBESGUEITERATOR_H

#include "MRGrid.h"
#include "GridNode.h"
#include "TreeIterator.h"

template<int D>
class LebesgueIterator: public TreeIterator<D> {
public:
    LebesgueIterator(MRGrid<D> *tree, int dir = TreeIterator<D>::TopDown):
        TreeIterator<D>(tree, dir) {
        init(tree);
    }
    virtual ~LebesgueIterator() {}
    inline bool next() {
        if (not this->state) {
	    return false;
	}
        if (this->mode == TreeIterator<D>::TopDown) {
            if (this->tryNode()) {
		return true;
	    }
	}
	GridNode<D> &node = *this->state->node;
	if ((node.getDepth() < this->maxDepth) and
		not (node.isEndNode() and not this->returnGenNodes)) {
	    const int nChildren = 1 << D;
	    for (int i = 0; i < nChildren; i++) {
		if (this->tryChild(i)) {
		    return true;
		}
	    }
	}
        if (this->tryNextRoot()) {
	    return true;
	}
        if (this->mode == TreeIterator<D>::BottomUp) {
            if (this->tryNode()) {
	    	return true;
	    }
	}
        this->removeState();
        return next();
    }
};

#endif // LEBESGUEITERATOR_H
