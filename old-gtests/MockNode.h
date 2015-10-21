#ifndef MOCKNODE_H
#define MOCKNODE_H

#include "MRNode.h"
#include "MockTree.h"

template<int D> class MockTree;

template<int D>
class MockNode: public MRNode<D> {
public:
    MockNode() : MRNode<D>() { }
    MockNode(MRTree<D> *t, const NodeIndex<D> &idx) : MRNode<D>(t, idx) { }
    MockNode(MRNode<D> *p, const int *l = 0) : MRNode<D>(p, l) { }
    MockNode(const MRNode<D> &nd, MRNode<D> *p, bool cc = true) : MRNode<D>(nd, p, cc) { }
    MockNode(const MRNode<D> &nd, MRTree<D> *t, bool cc = true) : MRNode<D>(nd, t, cc) { }
    virtual ~MockNode() { }

    MockNode<D> &operator=(const MRNode<D> &nd) {
        NOT_IMPLEMENTED_ABORT;
    }
protected:
    MRNode<D> *retrieveNode(int n, const double *r) {
        NOT_REACHED_ABORT;
    }
    MRNode<D> *retrieveNode(const NodeIndex<D> &idx) {
        NOT_REACHED_ABORT;
    }
    void createChild(int i) {
        NOT_REACHED_ABORT;
    }
    void genChild(int i) {
        NOT_REACHED_ABORT;
    }
};


#endif // MOCKNODE_H
