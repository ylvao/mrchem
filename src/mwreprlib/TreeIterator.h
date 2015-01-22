#ifndef TREEITERATOR_H
#define TREEITERATOR_H

#include "MRTree.h"
#include "MRNode.h"
#include "constants.h"

template <int D> struct IteratorNode;

template<int D>
class TreeIterator {
public:
    TreeIterator(MRTree<D> *tree, int dir = TopDown) {
        this->mode = dir;
        this->state = 0;
        this->initialState = 0;
        this->returnGenNodes = true;
    }
    virtual ~TreeIterator() {
        if (this->initialState != 0) {
            delete this->initialState;
        }
    }

    virtual bool next() = 0;

    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }
    void setMaxDepth(int depth) {
        if (depth < 0) {
            MSG_ERROR("Cannot have negative depth");
        }
        this->maxDepth = depth;
    }

    inline MRNode<D> &getNode() { return *this->state->node; }
    inline MRNode<D> &operator()() { return *this->state->node; }

    friend class IteratorNode<D>;
protected:
    int root;
    int nRoots;
    int mode;
    int maxDepth;
    bool returnGenNodes;
    IteratorNode<D> *state;
    IteratorNode<D> *initialState;

    void init(MRTree<D> *tree) {
        this->root = 0;
        this->maxDepth = tree->getMaxDepth();
        this->nRoots = tree->getRootBox().getNBoxes();
        this->state = new IteratorNode<D>(&tree->getRootBox().getNode(root));
        // Save the first state so it can be properly deleted later
        this->initialState = this->state;
    }
    bool tryNode() {
        if (not this->state) {
            return false;
        }
        if (this->state->doneNode) {
            return false;
        }
        this->state->doneNode = true;
        return true;
    }
    bool tryChild(int i) {
        if (not this->state) {
            return false;
        }
        if (this->state->doneChild[i]) {
            return false;
        }
        this->state->doneChild[i] = true;
        if (this->state->node->isLeafNode()) {
            return false;
        }
        if (this->state->node->isLeafNode()) {
            return false;
        }
        MRNode<D> *child = &this->state->node->getChild(i);
        this->state = new IteratorNode<D>(child, this->state);
        return next();
    }
    bool tryNextRoot() {
        if (not this->state) {
            return false;
        }
        if (not this->state->node->isRoot()) {
            return false;
        }
        this->root++;
        if (this->root >= this->nRoots) {
            return false;
        }
        MRNode<D> *nextRoot = &state->node->getTree().getRootBox().getNode(root);
        this->state = new IteratorNode<D>(nextRoot, this->state);
        return next();
    }
    void removeState() {
        if (this->state == this->initialState) {
            this->initialState = 0;
        }
        if (this->state != 0) {
            IteratorNode<D> *spare = this->state;
            this->state = spare->next;
            spare->next = 0;
            delete spare;
        }
    }
    void setDirection(int dir) {
        switch(dir) {
        case TopDown:
            this->mode = TopDown;
            break;
        case BottomUp:
            this->mode = BottomUp;
            break;
        default:
            MSG_FATAL("Invalid recursive direction!");
            break;
        }
    }
};

template<int D>
class IteratorNode {
public:
    MRNode<D> *node;
    IteratorNode<D> *next;
    bool doneNode;
    bool doneChild[1 << D];

    IteratorNode(MRNode<D> *_node, IteratorNode<D> *_next = 0) : node(_node), next(_next) {
        this->doneNode = false;
        int nChildren = 1 << D;
        for (int i = 0; i < nChildren; i++) {
            this->doneChild[i] = false;
        }
    }

    ~IteratorNode() { delete this->next; }
};

#endif // TREEITERATOR_H
