/**
 *  Simple n-dimensional node
 *
 *  Created on: May 29, 2009
 *      Author: jonas
 */

#ifndef MRNODE_H_
#define MRNODE_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>

#include "TelePrompter.h"
#include "macros.h"
#include "parallel.h"
#include "mwrepr_declarations.h"
#include "MRTree.h"
#include "HilbertPath.h"

#ifdef OPENMP
#define SET_NODE_LOCK() omp_set_lock(&this->node_lock)
#define UNSET_NODE_LOCK() omp_unset_lock(&this->node_lock)
#define TEST_NODE_LOCK() omp_test_lock(&this->node_lock)
#else
#define SET_NODE_LOCK()
#define UNSET_NODE_LOCK()
#define TEST_NODE_LOCK() false
#endif

template<int D>
class MRNode {
public:
    MRNode(MRTree<D> &t, const NodeIndex<D> &nIdx);
    MRNode(MRNode<D> &p, int cIdx);
    MRNode<D> &operator=(const MRNode<D> &nd);
    virtual ~MRNode();

    int getTDim() const { return this->tree->getTDim(); }
    int getKp1() const { return this->tree->getKp1(); }
    int getKp1_d() const { return this->tree->getKp1_d(); }
    int getOrder() const { return this->tree->getOrder(); }
    int getDepth() const { return this->nodeIndex->getScale()-this->tree->getRootScale(); }
    int getScale() const { return this->nodeIndex->getScale(); }
    int getRankId() const { return this->nodeIndex->getRankId(); }
    int getNChildren() const { if (isBranchNode()) return getTDim(); return 0; }

    const MRTree<D> &getMRTree() const {return *this->tree; }
    const MRNode<D> &getMRParent() const { return *this->parent; }
    const MRNode<D> &getMRChild(int i) const {
        assert(this->children != 0);
        return *this->children[i];
    }

    MRTree<D> &getMRTree() { return *this->tree; }
    MRNode<D> &getMRParent() { return *this->parent; }
    MRNode<D> &getMRChild(int i) {
        assert(this->children != 0);
        return *this->children[i];
    }

    const int *getTranslation() const { return nodeIndex->getTranslation(); }
    const NodeIndex<D> &getNodeIndex() const { return *nodeIndex; }

    void getCenter(double *r) const;
    void getLowerBounds(double *r) const;
    void getUpperBounds(double *r) const;

    inline bool isRoot() const;
    inline bool hasCoefs() const;
    inline bool isEndNode() const;
    inline bool isGenNode() const;
    inline bool isLeafNode() const;
    inline bool isAllocated() const;
    inline bool isBranchNode() const;
    inline bool hasChild(int i) const;

    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }
    void setHasCoefs(bool flag = true) {
        if (flag) {
            SET_BITS(status, FlagHasCoefs | FlagAllocated);
        } else {
            CLEAR_BITS(status, FlagHasCoefs);
        }
    }

    bool hasCoord(const double *r) const;
    bool isCompatible(const MRNode<D> &node);
    bool isAncestor(const NodeIndex<D> &idx) const;
    bool isDecendant(const NodeIndex<D> &idx) const;

    void calcChildTranslation(int *transl, int cIdx) const;
    void calcChildNodeIndex(NodeIndex<D> &nIdx, int cIdx) const;
    int getChildIndex(const NodeIndex<D> &nIdx) const;
    int getChildIndex(const double *r) const;

    virtual void genChildren();
    virtual void createChildren();
    virtual void deleteChildren();

    void lockNode() { SET_NODE_LOCK(); }
    void unlockNode() { UNSET_NODE_LOCK(); }
    bool testLock() { return TEST_NODE_LOCK(); }

    void setRankId(int n) { this->nodeIndex->setRankId(n); }
    bool isLocal() const {
        if (this->getRankId() == this->tree->getRankId()) {
            return true;
        }
        return false;
    }
    bool isCommon() const {
        if (this->getRankId() < 0) {
            return true;
        }
        return false;
    }
    bool isForeign() const {
        if (isLocal() or isCommon()) {
            return false;
        }
        return true;
    }

    const HilbertPath<D> &getHilbertPath() const { return *this->hilbertPath; }

    friend class GridNode<D>;

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const MRNode<T> &nd);

protected:
    MRTree<D> *tree;
    MRNode<D> *parent;	    ///< Parent node
    MRNode<D> **children;   ///< 2^D children

    NodeIndex<D> *nodeIndex;
    HilbertPath<D> *hilbertPath;

    bool diffBranch(const MRNode<D> &rhs) const;
    inline bool checkStatus(unsigned char mask) const;

    virtual MRNode<D> *retrieveNode(int n, const double *r);
    virtual MRNode<D> *retrieveNode(const NodeIndex<D> &idx);

    const MRNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx) const;
    MRNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx);

    const MRNode<D> *retrieveNodeOrEndNode(const double *r, int depth) const;
    const MRNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx) const;
    MRNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx);
    MRNode<D> *retrieveNodeOrEndNode(const double *r, int depth);

    void allocKindergarten();
    virtual void createChild(int i) = 0;
    virtual void genChild(int i) = 0;

    void purgeGenerated();

    void assignDecendantTags(int rank);
    void broadcastCoefs(int src, mpi::communicator *comm  = 0);
    virtual mpi::request isendCoefs(int who, int tag, int comp = -1) = 0;
    virtual mpi::request ireceiveCoefs(int who, int tag, int comp = -1) = 0;

    static const unsigned char FlagBranchNode =	B8(00000001);
    static const unsigned char FlagGenNode =	B8(00000010);
    static const unsigned char FlagHasCoefs =	B8(00000100);
    static const unsigned char FlagAllocated =	B8(00001000);
    static const unsigned char FlagEndNode =	B8(00010000);
    static const unsigned char FlagRootNode =	B8(00100000);
#ifdef OPENMP
    omp_lock_t node_lock;
#endif

    friend class MRTree<D>;
private:
    unsigned char status;

    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
        NOT_IMPLEMENTED_ABORT
                /*
          ar & parent;
          ar & tree;
          ar & nodeIndex;
          ar & hilbertPath;
          ar & status;
          if (this->isBranchNode()) {
           assert(children != 0);
           for (int i=0; i < tDim; i++) {
            ar & children[i];
           }
          }
          if (this->isEndNode() and not this->isForeign()) {
           if (checkStatus(this->FlagAllocated | this->FlagHasCoefs)) {
            const double *data = this->getRawData();
            for (int i = 0; i < this->getNCoefs(); i++) {
             ar & data[i];
            }
           }
          }
        */
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT
                /*
          ar & parent;
          ar & tree;
          ar & nodeIndex;
          ar & hilbertPath;
          ar & status;
          if (this->isBranchNode()) {
           allocKindergarten();
           for (int i=0; i < tDim; i++) {
            ar & children[i];
           }
          }
          if (this->isEndNode() and not this->isForeign()) {
           if (this->hasCoefs()) {
            this->allocCoefs();
            this->setHasCoefs(); // allocCoefs() resets hasCoefs!
            double *data = this->getCoefs().data();
            for (int i = 0; i < this->getNCoefs(); i++) {
             ar & data[i];
            }
            this->calcNorms();
           }
          } else {
           this->clearIsAllocated();
           this->setHasCoefs(false);
          }
          clearNodeWeight(0);
          clearNodeWeight(1);
        */
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

/** Allocation status of s/d-coefs is stored in the status bits for
 * serialization purposes. It's not enough to test if coefs == 0.
 */
template<int D>
bool MRNode<D>::isAllocated() const {
    if (this->status & FlagAllocated) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::hasCoefs() const {
    if (this->status & FlagHasCoefs) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::isGenNode() const {
    if (this->status & FlagGenNode) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::isLeafNode() const {
    if (this->status & FlagBranchNode) {
        return false;
    }
    return true;
}

template<int D>
bool MRNode<D>::isBranchNode() const {
    if (this->status & FlagBranchNode) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::isEndNode() const {
    if (this->status & FlagEndNode) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::isRoot() const {
    if (this->status & FlagRootNode) {
        return true;
    }
    return false;
}

template<int D>
bool MRNode<D>::hasChild(int i) const {
    assert(i >= 0 and i < getTDim());
    assert(this->children != 0);
    if (this->children[i] == 0) {
        return false;
    }
    return true;
}

template<int D>
bool MRNode<D>::checkStatus(unsigned char mask) const {
    if (mask == (this->status & mask)) {
        return true;
    }
    return false;
}

template<int D>
std::ostream& operator<<(std::ostream &o, const MRNode<D> &nd) {
    std::string flags ="      ";
    o << nd.nodeIndex;
    if (nd.isRoot()) {
        flags[0] = 'R';
    }
    if (nd.isEndNode()) {
        flags[1] = 'E';
    }
    if (nd.isBranchNode()) {
        flags[2] = 'B';
    } else {
        flags[2] = 'L';
    }
    if (nd.isGenNode()) {
        flags[3] = 'G';
    } else {
        flags[3] = 'P';
    }
    if (nd.isAllocated()) {
        flags[4] = 'A';
    }
    if (nd.hasCoefs()) {
        flags[5] = 'C';
    }
    o << " " << flags;
    return o;
}

#endif /* MRNODE_H_ */
