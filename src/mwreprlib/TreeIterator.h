#ifndef TREEITERATOR_H
#define TREEITERATOR_H

#include "MRGrid.h"
#include "GridNode.h"

template <int D> struct IteratorNode;

template<int D>
class TreeIterator {
public:
	TreeIterator(MRGrid<D> *tree, int dir = TopDown) {
		mode = dir;
		state = 0;
		initialState = 0;
		returnGenNodes = true;
	}
	virtual ~TreeIterator() {
		if (initialState != 0) {
			delete initialState;
		}
	}

	virtual bool next() = 0;

	void setReturnGenNodes(bool i = true) {
		returnGenNodes = i;
	}
	void setMaxDepth(int depth) {
		if (depth < 0) {
			MSG_ERROR("Cannot have negative depth");
		}
		maxDepth = depth;
	}
	inline GridNode<D> &getNode() {
		return *state->node;
	}
	inline GridNode<D> &operator()() {
		return *state->node;
	}

	enum {
		TopDown,
		BottomUp
	};

	friend class IteratorNode<D>;
protected:
	int root;
	int nRoots;
	int mode;
	int maxDepth;
	bool returnGenNodes;
	IteratorNode<D> *state;
	IteratorNode<D> *initialState;

	void init(MRGrid<D> *tree) {
		root = 0;
		maxDepth = tree->getMaxDepth();
		nRoots = tree->getRootBox().getNBoxes();
		state = new IteratorNode<D>(&tree->getRootBox().getNode(root));
		// Save the first state so it can be properly deleted later
		initialState = state;
	}
	bool tryNode() {
		if (not state) {
			return false;
		}
		if (state->doneNode) {
			return false;
		}
		state->doneNode = true;
		return true;
	}
	bool tryChild(int i) {
		if (not state) {
			return false;
		}
		if (state->doneChild[i]) {
			return false;
		}
		state->doneChild[i] = true;
		if (state->node->isLeafNode()) {
			return false;
		}
		if (state->node->isLeafNode()) {
			return false;
		}
		state = new IteratorNode<D>(state->node->getChild(i), state);
		return next();
	}
	bool tryNextRoot() {
		if (not state) {
			return false;
		}
		if (not state->node->isRoot()) {
			return false;
		}
		root++;
		if (root >= nRoots) {
			return false;
		}
		state = new IteratorNode<D>(
				&state->node->getGrid().getRootBox().getNode(root), state);
		return next();
	}
	void removeState() {
		if (state == initialState) {
			initialState = 0;
		}
		if (state != 0) {
			IteratorNode<D> *spare = state;
			state = spare->next;
			spare->next = 0;
			delete spare;
		}
	}
	void setDirection(int dir) {
		switch(dir) {
		case TopDown:
			mode = TopDown;
			break;
		case BottomUp:
			mode = BottomUp;
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
	GridNode<D> *node;
	IteratorNode<D> *next;
	bool doneNode;
	bool doneChild[1 << D];

	IteratorNode(GridNode<D> *_node, IteratorNode<D> *_next = 0)
		: node(_node), next(_next) {
		doneNode = false;
		int nChildren = 1 << D;
		for (int i = 0; i < nChildren; i++) {
			doneChild[i] = false;
		}
	}

	~IteratorNode() { delete next; }
};

#endif // TREEITERATOR_H
