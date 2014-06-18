#include "MRTree.h"
#include "MathUtils.h"
#include "MRNode.h"
#include "LebesgueIterator.h"


using namespace std;
using namespace Eigen;

template<int D> int MRTree<D>::defaultMaxDepth = 30;
template<int D> int MRTree<D>::defaultOrder = 3;

template<int D>
MRTree<D>::MRTree(int k, const NodeBox<D> *box) {
    if (k < 0) {
	k = defaultOrder;
    }
    this->order = k;
    this->kp1 = this->order + 1;
    this->kp1_d = MathUtils::ipow(this->kp1, D);

    if (box != 0) {
        this->rootBox = new NodeBox<D>(*box);
    } else {
	NOT_IMPLEMENTED_ABORT
    }

    this->maxDepth = defaultMaxDepth;
    this->maxScale = this->getRootScale() + this->maxDepth - 1;

    this->nNodes = 0;
    this->nodesAtDepth.push_back(0);

    mpi::communicator world;
    this->rank = world.rank();
    this->nThreads = omp_get_max_threads();
    allocNodeCounters();

#ifdef OPENMP
    omp_init_lock(&tree_lock);
#endif
}

template<int D>
MRTree<D>::MRTree(const MRTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
MRTree<D>::~MRTree() {
    this->endNodeTable.clear();
    delete this->rootBox;
    if (this->nNodes != 0) {
	THROW_ERROR("Node count != 0 -> " << this->nNodes);
    }
    if (this->nodesAtDepth.size() != 0) {
	THROW_ERROR("Nodes at depth != 0 -> " << this->nodesAtDepth.size());
    }
    deleteNodeCounters();
#ifdef OPENMP
    omp_destroy_lock(&tree_lock);
#endif
}

template<int D>
void MRTree<D>::allocNodeCounters() {
    this->nGenNodes = new int[this->nThreads];
    this->nAllocGenNodes = new int[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->nGenNodes[i] = 0;
        this->nAllocGenNodes[i] = 0;
    }
}

template<int D>
void MRTree<D>::deleteNodeCounters() {
    delete[] this->nGenNodes;
    delete[] this->nAllocGenNodes;
}

template<int D>
void MRTree<D>::copyTreeParams(const MRTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRTree<D>::initializeRootNodes() {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRTree<D>::yieldChildren(MRNodeVector &nodeTable, const NodeIndexSet &idxSet) {
    typename set<const NodeIndex<D> *>::iterator it;
    for (it = idxSet.begin(); it != idxSet.end(); it++) {
        MRNode<D> &node = getNode(**it);
        int childDepth = node.getDepth() + 1;
        if (this->maxDepth != 0 and childDepth > this->maxDepth) {
            println(1, "+++ Maximum depth reached: " << childDepth);
            node.setIsEndNode();
        } else {
            node.createChildren();
            for (int i = 0; i < node.getNChildren(); i++) {
                MRNode<D> *child = &node.getChild(i);
                nodeTable.push_back(child);
            }
        }
    }
}

template<int D>
void MRTree<D>::setDefaultOrder(int k) {
    assert(k > 0);
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRTree<D>::setDefaultMaxDepth(int max_depth) {
    assert(max_depth > 0);
    NOT_IMPLEMENTED_ABORT
}

template<int D>
bool MRTree<D>::diffTree(const MRTree<D> &tree) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
bool MRTree<D>::checkCompatible(const MRTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MRTree<D>::incrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    assert(depth <= this->maxDepth);
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
        for (int i = 0; i < depth - n; i++) {
            this->nodesAtDepth.push_back(0);
        }
    }
    this->nodesAtDepth[depth]++;
    this->nNodes++;
}

template<int D>
void MRTree<D>::decrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    assert(depth < this->nodesAtDepth.size());
    this->nodesAtDepth[depth]--;
    assert(this->nodesAtDepth[depth] >= 0);
    if (this->nodesAtDepth[depth] == 0) { 
        this->nodesAtDepth.pop_back();
    }

    this->nNodes--;
    assert(this->nNodes >= 0);
}

template<int D>
void MRTree<D>::updateGenNodeCounts() {
    lockTree();
    for (int i = 1; i < this->nThreads; i++) {
        this->nGenNodes[0] += this->nGenNodes[i];
        this->nAllocGenNodes[0] += this->nAllocGenNodes[i];
        this->nGenNodes[i] = 0;
        this->nAllocGenNodes[i] = 0;
    }
    assert(this->nGenNodes[0] >= 0);
    assert(this->nAllocGenNodes[0] >= 0);
    unlockTree();
}

template<int D>
void MRTree<D>::incrementGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]++;
}

template<int D>
void MRTree<D>::decrementGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]--;
}

template<int D>
void MRTree<D>::incrementAllocGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nAllocGenNodes[n]++;
}

template<int D>
void MRTree<D>::decrementAllocGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nAllocGenNodes[n]--;
}

template<int D>
int MRTree<D>::getNNodes(int depth) const {
    if (depth < 0) {
	return this->nNodes;
    }
    if (depth >= this->nodesAtDepth.size()) {
	return 0;
    }
    return this->nodesAtDepth[depth];
}

template<int D>
int MRTree<D>::getNAllocGenNodes() {
    updateGenNodeCounts();
    return this->nAllocGenNodes[0];
}

template<int D>
int MRTree<D>::getNGenNodes() {
    updateGenNodeCounts();
    return this->nGenNodes[0];
}

template<int D>
const MRNode<D>* MRTree<D>::findNode(const NodeIndex<D> &idx) const {
    assert(idx.getScale() <= getMaxScale());
    const MRNode<D> &root = getRootNode(idx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

template<int D>
MRNode<D>* MRTree<D>::findNode(const NodeIndex<D> &idx) {
    assert(idx.getScale() <= getMaxScale());
    MRNode<D> &root = getRootNode(idx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

template<int D>
const MRNode<D>* MRTree<D>::findNode(const double *r, int depth) const {
    const MRNode<D> &root = getRootNode(r);
    println(0, "Root " << root);
    return root.retrieveNodeOrEndNode(r, depth);
}

template<int D>
MRNode<D>* MRTree<D>::findNode(const double *r, int depth) {
    MRNode<D> &root = getRootNode(r);
    return root.retrieveNodeOrEndNode(r, depth);
}

template<int D>
MRNode<D>& MRTree<D>::getNode(const NodeIndex<D> &idx) {
    assert(idx.getScale() <= getMaxScale());
    MRNode<D> &root = getRootNode(idx);
    assert(root.isAncestor(idx));
    return *(root.retrieveNode(idx));
}

template<int D>
MRNode<D>& MRTree<D>::getNode(const double *r, int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>& MRTree<D>::getNodeNoGen(const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::makeNodeTable(MRNodeVector &nodeTable) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::makeNodeTable(std::vector<MRNodeVector > &nodeTable) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::makeLocalNodeTable(MRNodeVector &nodeTable) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::makeLocalNodeTable(std::vector<MRNodeVector > &nodeTable) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::copyEndNodeTable(MRNodeVector &nodeTable) {
    for (int n = 0; n < getNEndNodes(); n++) {
	MRNode<D> &node = getEndNode(n);
	nodeTable.push_back(&node);
    }
}

template<int D>
void MRTree<D>::resetEndNodeTable() {
    clearEndNodeTable();
    LebesgueIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MRNode<D> &node = it.getNode();
        if (node.isEndNode()) {
	    this->endNodeTable.push_back(&node);
        }
    }
}

template<int D>
void MRTree<D>::purgeGenerated() {
    for (int n = 0; n < getNEndNodes(); n++) {
        getEndNode(n).purgeGenerated();
    }
} 

template<int D>
int MRTree<D>::countBranchNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
int MRTree<D>::countLeafNodes(int depth) {
    int nNodes = 0;
    LebesgueIterator<D> it(this);
    while (it.next()) {
        MRNode<D> &node = it.getNode();
	if (node.getDepth() == depth or depth < 0) {
            if (node.isLeafNode()) {
		nNodes++;
	    }
	}
    }
    return nNodes;
}

template<int D>
int MRTree<D>::countAllocNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
int MRTree<D>::countMyNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::printNodeRankCount() {
    mpi::communicator world;
    int mDepth = getDepth();
    int count[mDepth + 1][world.size()];
    for (int i = 0; i < mDepth + 1; i++) {
        for (int j = 0; j < world.size(); j++) {
            count[i][j] = 0;
        }
    }
    LebesgueIterator<D> it(this);
    int notDistributed = 0;
    while(it.next()) {
        MRNode<D> &node = it.getNode();
        if (node.isGenNode()) {
            continue;
        }
        int depth = node.getDepth();
        int rank = node.getRankId();
	if (rank >= 0) {
            count[depth][rank]++;
	    count[mDepth][rank]++;
	} else {
	    notDistributed++;
	}
    }

    println(0, endl);
    printout(0, "   |");
    for (int j = 0; j < world.size(); j++) {
        printout(0, setw(5) << j);
    }
    println(0, "|" << endl);
    for (int i = 0; i < mDepth + 1; i++) {
        if (i == mDepth) {
            printout(0, endl);
        }
        printout(0, setw(3) << i << "|");
        for (int j = 0; j < world.size(); j++) {
            printout(0, setw(5) << count[i][j]);
        }
        println(0, "|");
    }
    println(0, endl);
    if (notDistributed > 0) {
        println(0, "There were " << notDistributed << " non-distributed nodes");
	println(0, endl);
    }
}

template<int D>
void MRTree<D>::broadcastTree() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRTree<D>::distributeEndNodes() {
    mpi::communicator world;
    int nLocales = world.size();
    int nNodes = getNEndNodes();
    int tDim = getTDim();
    int nParents = nNodes/tDim;

    int *distNodes = new int[nLocales];
    for (int i = 0; i < nLocales; i++) {
	distNodes[i] = 0;
    }
    
    int rank = 0;
    int nPerRank = nNodes/nLocales;
    println(0, "Number of nodes       " << nNodes);
    println(0, "Parents to distribute " << nParents);
    println(0, "Number of MPI hosts   " << nLocales);
    println(0, "Nodes per host        " << nPerRank);
    for (int n = 0; n < nNodes; n++) {
	MRNode<D> &parent = getEndNode(n).getParent();
	for (int i = 0; i < tDim; i++) {
	    MRNode<D> &child = parent.getChild(i);
	    if (child.isEndNode() and child.getRankId() < 0) {
		assert(rank < nLocales);
	        child.setRankId(rank);
	        distNodes[rank]++;
	    }
	}
	if (distNodes[rank] >= nPerRank) {
	    rank++;
	}
    }
    int nDist = 0;
    for (int i = 0; i < nLocales; i++) {
	println(0, "MPI rank " << i << " has " << distNodes[i] << " nodes");
	nDist += distNodes[i];
    }
    if (nDist != nNodes) {
	THROW_ERROR("Not all endNodes were distributed");
    }
}

template class MRTree<1>;
template class MRTree<2>;
template class MRTree<3>;

/*
 
template<int D>
void MRGrid<D>::incrementNodeCount(int scale) {
    int depth = scale - this->getRootScale();
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
        for (int i = 0; i < depth - n; i++) {
            this->nodesAtDepth.push_back(0);
        }
    }
    int nodes = this->nodesAtDepth[depth];
    nodes++;
    this->nodesAtDepth[depth] = nodes;
}

template<int D>
void MRGrid<D>::decrementNodeCount(int scale) {
    unsigned int depth = scale - this->getRootScale();
    if (depth >= this->nodesAtDepth.size()) {
        THROW_ERROR("Depth out of bounds!");
    }
    int nodes = this->nodesAtDepth[depth];
    nodes--;
    if (nodes < 0) {
        THROW_ERROR("Number of nodes cannot be negative.");
    }
    this->nodesAtDepth[depth] = nodes;
    if (nodes == 0) { 
        this->nodesAtDepth.pop_back();
    }
}

template<int D>
int MRGrid<D>::getNodesAtDepth(int depth) const {
    if (depth < 0) {
	THROW_ERROR("Negative depth");
    }
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
	return 0;
    }
    return this->nodesAtDepth[depth];
}

template<int D>
void MRGrid<D>::clearEndNodeTable() {
    this->endNodeTable.clear();
}

template<int D>
void MRGrid<D>::resetEndNodeTable() {
    this->endNodeTable.clear();
    LebesgueIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
        if (node.isEndNode()) {
	    this->endNodeTable.push_back(&node);
        }
    }
}

template<int D>
void MRGrid<D>::copyEndNodeTable(GridNodeVector &outTable) {
    for (int i = 0; i < this->endNodeTable.size(); i++) {
	GridNode<D> *node = this->endNodeTable[i];
	outTable.push_back(node);
    }
}
    
template<int D>
int MRGrid<D>::getNNodes(int depth) const {
    int nScales = this->nodesAtDepth.size();
    if (depth > nScales) {
	return 0;
    }
    if (depth >= 0) {
	return this->nodesAtDepth[depth];
    }
    int totNodes = 0;
    for (int i = 0; i < nScales; i++) {
	totNodes += this->nodesAtDepth[i];
    }
    return totNodes; 
}

template<int D>
int MRGrid<D>::countBranchNodes(int depth) {
    int nNodes = 0;
    LebesgueIterator<D> it(this);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
	if (node.getDepth() == depth or depth < 0) {
            if (node.isBranchNode()) {
		nNodes++;
	    }
	}
    }
    return nNodes;
}

template<int D>
int MRGrid<D>::countLeafNodes(int depth) {
    int nNodes = 0;
    LebesgueIterator<D> it(this);
    while (it.next()) {
        GridNode<D> &node = it.getNode();
	if (node.getDepth() == depth or depth < 0) {
            if (node.isLeafNode()) {
		nNodes++;
	    }
	}
    }
    return nNodes;
}

template<int D>
void MRGrid<D>::yieldChildren(GridNodeVector &nodeTable, const NodeIndexSet &idxSet) {
    typename set<const NodeIndex<D> *>::iterator it;
    for (it = idxSet.begin(); it != idxSet.end(); it++) {
        GridNode<D> &node = this->rootBox->getNode(**it);
        int childDepth = node.getDepth() + 1;
        if (this->maxDepth != 0 and childDepth > this->maxDepth) {
            println(1, "+++ Maximum depth reached: " << childDepth);
            node.setIsEndNode();
        } else {
            node.createChildren();
            for (int i = 0; i < node.getNChildren(); i++) {
                GridNode<D> *child = node.getChild(i);
                nodeTable.push_back(child);
            }
        }
    }
}
*/
