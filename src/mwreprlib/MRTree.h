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

#include <set>
#include <vector>

#include "NodeBox.h"
#include "parallel.h"
#include "mwrepr_declarations.h"

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
    MRTree(int k = defaultOrder, const NodeBox<D> *box = 0);
    MRTree(const MRTree<D> &tree);
    virtual ~MRTree();

    virtual void clear() = 0;
    void copyTreeParams(const MRTree<D> &tree);

    void setName(const std::string &n) { this->name = n; }
    const std::string &getName() const { return this->name; }

    int getKp1() const { return this->kp1; }
    int getKp1_d() const { return this->kp1_d; }
    int getOrder() const { return this->order; }
    int getDim() const { return D; }
    int getTDim() const { return this->tDim; }
    int getNNodes(int depth = -1) const;
    int getMaxDepth() const { return this->maxDepth; }
    int getMaxScale() const { return this->maxScale; }
    int getNEndNodes() const { return this->endNodeTable.size(); }
    int getNAllocGenNodes();
    int getNGenNodes();
    virtual int getDepth() const { return this->nodesAtDepth.size(); }

    NodeBox<D> &getRootBox() { return *this->rootBox; }
    const NodeBox<D> &getRootBox() const { return *this->rootBox; }
    int getRootScale() const { return this->rootBox->getRootScale(); }
    int getNRootNodes() const { return this->rootBox->getNBoxes(); }
    const double *getOrigin() const { return this->rootBox->getOrigin(); }
    const double *getLowerBounds() const { return this->rootBox->getLowerBounds(); }
    const double *getUpperBounds() const { return this->rootBox->getUpperBounds(); }

    const MRNode<D> *findNode(const double *r, int depth = -1) const;
    const MRNode<D> *findNode(const NodeIndex<D> &idx) const;
    MRNode<D> *findNode(const NodeIndex<D> &idx);
    MRNode<D> *findNode(const double *r, int depth = -1);
    MRNode<D> &getNode(const NodeIndex<D> &idx);
    MRNode<D> &getNode(const double *r, int depth = -1);
    MRNode<D> &getNodeNoGen(const NodeIndex<D> &idx);

    const MRNode<D> &getEndNode(int i) const { return *this->endNodeTable[i]; }
    const MRNode<D> &getRootNode(int i = 0) const { return this->rootBox->getNode(i); }
    const MRNode<D> &getRootNode(const double *r) const { return this->rootBox->getNode(r); }
    const MRNode<D> &getRootNode(const NodeIndex<D> &idx) const { return this->rootBox->getNode(idx); }

    MRNode<D> &getEndNode(int i) { return *this->endNodeTable[i]; }
    MRNode<D> &getRootNode(int i = 0) { return this->rootBox->getNode(i); }
    MRNode<D> &getRootNode(const double *r) { return this->rootBox->getNode(r); }
    MRNode<D> &getRootNode(const NodeIndex<D> &idx) { return this->rootBox->getNode(idx); }

    void purgeGenerated();

    void broadcastTree();
    void distributeEndNodes();

    void lockTree() { SET_TREE_LOCK(); }
    void unlockTree() { UNSET_TREE_LOCK(); }
    bool testLock() { return TEST_TREE_LOCK(); }

    int getNThreads() const { return this->nThreads; }
    int getRankId() const { return this->rank; }
    bool isScattered() const { return this->scattered; }

    bool diffTree(const MRTree<D> &tree) const;
    bool checkCompatible(const MRTree<D> &tree);

    virtual bool saveTree(const std::string &file) = 0;
    virtual bool loadTree(const std::string &file) = 0;

    int countBranchNodes(int depth = -1);
    int countLeafNodes(int depth = -1);
    int countAllocNodes(int depth = -1);
    int countMyNodes(int depth = -1);
    void printNodeRankCount();

    static void setDefaultOrder(int _order);
    static void setDefaultMaxDepth(int max_depth);

    friend class MRNode<D>;
    friend class GridNode<D>;

protected:
    // Parameters that are set in construction and should never change
    int order;		
    int maxDepth;
    int rank;
    int nThreads;

    // Parameters that are derived internally and should not be set explicitly
    int kp1;
    int kp1_d;
    int maxScale;
    bool scattered;

    // Parameters that are dynamic and can be set by user
    std::string name;

    // Tree data
    int nNodes;
    int *nGenNodes;
    int *nAllocGenNodes;
    NodeBox<D> *rootBox;
    MRNodeVector endNodeTable;	    ///< Final projected node
    std::vector<int> nodesAtDepth;  ///< used to determine who has a given node

    const static int tDim = (1 << D);
    static int defaultSplitType;
    static int defaultScalingType;
    static int defaultUniformDepth;
    static int defaultMaxDepth;
    static int defaultOrder;

    int getRootIndex(const double *r) const { return this->rootBox->getBoxIndex(r); }
    int getRootIndex(const NodeIndex<D> &idx) { return this->rootBox->getBoxIndex(idx); }

    void allocNodeCounters();
    void deleteNodeCounters();

    void incrementNodeCount(int scale);
    void decrementNodeCount(int scale);
    void updateGenNodeCounts();
    void incrementGenNodeCount();
    void decrementGenNodeCount();
    void incrementAllocGenNodeCount();
    void decrementAllocGenNodeCount();

    virtual void initializeRootNodes() = 0;
    void yieldChildren(MRNodeVector &nodeTable, const NodeIndexSet &idxSet);

    void makeNodeTable(MRNodeVector &nodeTable);
    void makeNodeTable(std::vector<MRNodeVector > &nodeTable);

    void makeLocalNodeTable(MRNodeVector &nodeTable);
    void makeLocalNodeTable(std::vector<MRNodeVector > &nodeTable);

    void copyEndNodeTable(MRNodeVector &nodeTable);
    void resetEndNodeTable();
    void clearEndNodeTable() { this->endNodeTable.clear(); }

    void findMissingNodes(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);
    void findMissingParents(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);
    void findMissingChildren(MRNodeVector &nodeTable, std::set<MRNode<D> *> &missing);

    void tagNodes(MRNodeVector &nodeList, int _rank);

#ifdef OPENMP
    omp_lock_t tree_lock;
#endif
};

#endif /* MRTREE_H_ */
