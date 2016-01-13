/**
 *
 * \date Jun 5, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 */

#ifndef MRTREE_H_
#define MRTREE_H_

#include "parallel.h"
#include "mwrepr_declarations.h"

#include "NodeBox.h"

#ifdef OPENMP
#define SET_TREE_LOCK() omp_set_lock(&this->tree_lock)
#define UNSET_TREE_LOCK() omp_unset_lock(&this->tree_lock)
#define TEST_TREE_LOCK() omp_test_lock(&this->tree_lock)
#else
#define SET_TREE_LOCK()
#define UNSET_TREE_LOCK()
#define TEST_TREE_LOCK() false
#endif

template<int D>
class MRTree {
public:
    MRTree(const BoundingBox<D> &box);
    MRTree(const MRTree<D> &tree);
    virtual ~MRTree();

    void setName(const std::string &n) { this->name = n; }
    const std::string &getName() const { return this->name; }

    int getDim() const { return D; }
    int getTDim() const { return this->tDim; }
    int getNNodes(int depth = -1) const;
    int getNEndNodes() const { return this->endNodeTable.size(); }
    int getNAllocGenNodes();
    int getNGenNodes();
    int getRootScale() const { return this->rootBox.getScale(); }
    virtual int getDepth() const { return this->nodesAtDepth.size(); }

    NodeBox<D> &getRootBox() { return this->rootBox; }
    const NodeBox<D> &getRootBox() const { return this->rootBox; }

    MRNode<D> *findNode(const NodeIndex<D> &nIdx);
    const MRNode<D> *findNode(const NodeIndex<D> &nIdx) const;

    MRNode<D> &getNode(const NodeIndex<D> &nIdx);
    MRNode<D> &getNodeOrEndNode(const NodeIndex<D> &nIdx);
    const MRNode<D> &getNodeOrEndNode(const NodeIndex<D> &nIdx) const;

    MRNode<D> &getNode(const double *r, int depth = -1);
    MRNode<D> &getNodeOrEndNode(const double *r, int depth = -1);
    const MRNode<D> &getNodeOrEndNode(const double *r, int depth = -1) const;

    MRNode<D> &getEndNode(int i) { return *this->endNodeTable[i]; }
    const MRNode<D> &getEndNode(int i) const { return *this->endNodeTable[i]; }

    void lockTree() { SET_TREE_LOCK(); }
    void unlockTree() { UNSET_TREE_LOCK(); }
    bool testLock() { return TEST_TREE_LOCK(); }

    int getNThreads() const { return this->nThreads; }
    int getRankId() const { return this->rank; }

    virtual bool saveTree(const std::string &file) { NOT_IMPLEMENTED_ABORT; }
    virtual bool loadTree(const std::string &file) { NOT_IMPLEMENTED_ABORT; }

    int countBranchNodes(int depth = -1);
    int countLeafNodes(int depth = -1);
    int countAllocNodes(int depth = -1);
    int countMyNodes(int depth = -1);
    void printNodeRankCount();

    void checkGridOverlap(MRTree<D> &tree);
    void checkRankOverlap(MRTree<D> &tree);

    friend class MRNode<D>;
    friend class GenNode<D>;
protected:
    // Parameters that are set in construction and should never change
    const int rank;
    const int nThreads;

    // Parameters that are dynamic and can be set by user
    std::string name;

    // Tree data
    int nNodes;
    int *nGenNodes;
    int *nAllocGenNodes;
    NodeBox<D> rootBox;            ///< The actual container of nodes
    MRNodeVector endNodeTable;	   ///< Final projected nodes
    std::vector<int> nodesAtDepth;  ///< Node counter

    // Static default parameters
    const static int tDim = (1 << D);

    int getRootIndex(const double *r) const { return this->rootBox.getBoxIndex(r); }
    int getRootIndex(const NodeIndex<D> &nIdx) { return this->rootBox.getBoxIndex(nIdx); }

    void allocNodeCounters();
    void deleteNodeCounters();

    void incrementNodeCount(int scale);
    void decrementNodeCount(int scale);
    void updateGenNodeCounts();
    void incrementGenNodeCount();
    void decrementGenNodeCount();
    void incrementAllocGenNodeCount();
    void decrementAllocGenNodeCount();

    void splitNodes(const NodeIndexSet &idxSet, MRNodeVector *nVec = 0);

    void makeNodeTable(MRNodeVector &nodeTable);
    void makeNodeTable(std::vector<MRNodeVector > &nodeTable);

    void makeLocalNodeTable(MRNodeVector &nodeTable, bool common = false);
    void makeLocalNodeTable(std::vector<MRNodeVector > &nodeTable, bool common = false);

    MRNodeVector* copyEndNodeTable();
    MRNodeVector* getEndNodeTable() { return &this->endNodeTable; }

    void resetEndNodeTable();
    void clearEndNodeTable() { this->endNodeTable.clear(); }

    void findMissingNodes(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);
    void findMissingParents(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);
    void findMissingChildren(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);

    void distributeNodes(int depth = -1);
    void deleteForeign(bool keepEndNodes = false);

    void tagNodes(MRNodeVector &nodeList, int rank);
    void tagDecendants(MRNodeVector &nodeList);
    void distributeNodeTags(MRNodeVector &nodeList);

    void deleteGenerated();
    void clearGenerated();

#ifdef OPENMP
    omp_lock_t tree_lock;
#endif
};

#endif /* MRTREE_H_ */
