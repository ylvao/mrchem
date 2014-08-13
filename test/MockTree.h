#ifndef MOCKTREE_H
#define MOCKTREE_H

#include "MRTree.h"
#include "MockNode.h"

template<int D> class MockNode;

template<int D>
class MockTree: public MRTree<D> {
public:
    MockTree(int k, const BoundingBox<D> *box) : MRTree<D>(k, box) {
        initializeRootNodes();
        this->resetEndNodeTable();
    }
    MockTree(const MockTree<D> &tree) : MRTree<D>(tree) {
        this->rootBox = tree.rootBox;
        this->resetEndNodeTable();
    }
    virtual ~MockTree() { }

    virtual void clear() {
        NOT_IMPLEMENTED_ABORT;
    }
    virtual bool saveTree(const std::string &file) {
        NOT_IMPLEMENTED_ABORT;
        return false;
    }
    virtual bool loadTree(const std::string &file) {
        NOT_IMPLEMENTED_ABORT;
        return false;
    }
protected:
    virtual void initializeRootNodes() {
        for (int i = 0; i < this->getRootBox().getNBoxes(); i++) {
            NodeIndex<D> idx = this->getRootBox().getNodeIndex(i);
            MRNode<D> *root = new MockNode<D>(this, idx);
            this->rootBox->setNode(i, &root);
        }
    }
};


#endif // MOCKTREE_H
