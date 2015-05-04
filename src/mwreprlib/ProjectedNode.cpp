/**
 *
 *
 *  \date Aug 14, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 */

#include "ProjectedNode.h"
#include "GridNode.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

/** Root node constructor. By default root nodes are initialized to
  * represent functions which are constant zero. */
template<int D>
ProjectedNode<D>::ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx)
        : FunctionNode<D> (t, nIdx) {
    this->setIsEndNode();
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
        this->zeroNorms();
    }
}

/** ProjectedNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
ProjectedNode<D>::ProjectedNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx) {
    this->setIsEndNode();
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
        this->zeroNorms();
    }
}

/* Recurcive root node construcor*/
/*
template<int D>
ProjectedNode<D>::ProjectedNode(FunctionTree<D> &t, const GridNode<D> &gNode)
        : FunctionNode<D> (t, gNode.getNodeIndex()) {
    this->setIsEndNode();
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
        this->zeroNorms();
    }

    if (gNode.isBranchNode()) {
        this->allocKindergarten();
        this->setIsBranchNode();
        this->clearIsEndNode();
    }
    for (int cIdx = 0; cIdx < gNode.getNChildren(); cIdx++) {
        const GridNode<D> &gChild = gNode.getGridChild(cIdx);
        this->children[cIdx] = new ProjectedNode(*this, cIdx, gChild);
    }
}
*/
/* Recurcive node constructor*/
/*
template<int D>
ProjectedNode<D>::ProjectedNode(ProjectedNode<D> &p, int cIdx,
                                const MRNode<D> &gNode)
        : FunctionNode<D> (p, cIdx) {
    this->setIsEndNode();
    this->setRankId(gNode.getRankId());
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
        this->zeroNorms();
    }

    if (gNode.isBranchNode()) {
        this->allocKindergarten();
        this->setIsBranchNode();
        this->clearIsEndNode();
    }
    for (int cIdx = 0; cIdx < gNode.getNChildren(); cIdx++) {
        const GridNode<D> &gChild = gNode.getGridChild(cIdx);
        ProjectedNode<D> *child = new ProjectedNode(*this, cIdx, gChild);
        this->children[cIdx] = child;
    }
}
*/
/** ProjectedNode equals operator.
  * Copying the content of a node, not its location. Includes recursive copying
  * of children nodes. */
template<int D>
ProjectedNode<D> &ProjectedNode<D>::operator=(const ProjectedNode<D> &nd) {
    NOT_IMPLEMENTED_ABORT;
//    if (this == &node) {
//        return *this;
//    }
//    FunctionNode<D>::operator=(node);
//    copyRedundancyMap(node);
//    if (this->isLeafNode()) {
//        return *this;
//    }
//    for (int i = 0; i < this->tDim; i++) {
//        if (this->isEndNode()) {
//            const GenNode<D> &child =
//                    static_cast<const GenNode<D> &>(node.getMWChild(i));
//            this->children[i] = new GenNode<D>(child, this, true);
//        } else {
//            const ProjectedNode<D> &child =
//                    static_cast<const ProjectedNode<D> &>(node.getMWChild(i));
//            this->children[i] = new ProjectedNode<D>(child, this, true);
//        }
//    }
//    return *this;
}

/** Allocating children nodes.
  *
  * This routine creates 2^D empty ProjectedNode children nodes with the
  * appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChildren() {
    if (this->children == 0) {
        this->allocKindergarten();
    }
    for (int i = 0; i < this->getTDim(); i++) {
        createChild(i);
    }
    this->setIsBranchNode();
    this->clearIsEndNode();
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty ProjectedNode child node
  * with the appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChild(int i) {
    assert(this->children != 0);
    assert(this->children[i] == 0);
    ProjectedNode<D> *child = new ProjectedNode<D>(*this, i);
    this->children[i] = child;
}

/** Generating children nodes.
  *
  * This routine creates 2^D children GenNodes with the appropriate translation
  * and Hilbert path parameters. Option to calculate the children scaling coefs
  * by MW transform. */
template<int D>
void ProjectedNode<D>::genChildren(bool genEmpty) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->children == 0) {
//        this->allocKindergarten();
//    }
//    this->setIsBranchNode();
//    for (int i = 0; i < this->tDim; i++) {
//        genChild(i);
//    }
//    // On dist trees this is not always the case
//    if (this->hasCoefs() and not genEmpty) {
//        this->giveChildrenScaling();
//        for (int i = 0; i < this->tDim; i++) {
//            this->children[i]->calcNorms();
//        }
//    }
}

/** Generating child node.
  *
  * This routine creates 2^D children empty GenNodes with the appropriate
  * translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::genChild(int i) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->children[i] == 0);
//    int l[D];
//    this->calcChildTranslation(i, l);
//    FunctionNode<D> *child = new GenNode<D> (this, l);
//    this->children[i] = child;
//    return *child;
}

/** Given the scaling and wavelet coefficients of the node, do a wavelet
  * transform and distribute scaling coefficient to all children. Option to
  * overwrite or add up existing coefficients. */
//template<int D>
//void ProjectedNode<D>::giveChildrenScaling(bool overwrite) {
//    NOT_IMPLEMENTED_ABORT;
//    assert(this->isBranchNode());
//    if (not this->hasCoefs()) {
//        MSG_FATAL("No coefficients! You need to project first!")
//    }
//    ProjectedNode<D> tmpNode(*this, this->tree);
//    tmpNode.mwTransform(Reconstruction);
//    int kp1_d = this->getKp1_d();
//    for (int i = 0; i < this->tDim; i++) {
//        MWNode<D> &child = this->getMWChild(i);
//        if (not child.hasCoefs()) {
//            child.setCoefs(tmpNode.getCoefs().segment(i * kp1_d, kp1_d));
//        } else if (overwrite) {
//            child.getCoefs().segment(0, kp1_d) =
//                tmpNode.getCoefs().segment(i * kp1_d, kp1_d);
//        } else {
//            child.getCoefs().segment(0, kp1_d) +=
//                tmpNode.getCoefs().segment(i * kp1_d, kp1_d);
//        }
//        child.calcNorms();
//    }
//    tmpNode.clearTreePointer();
//}


/** Calculate all 2^D component norms (NOT squared norms!)*/
template<int D>
void ProjectedNode<D>::calcComponentNorms() {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->hasCoefs());
//    if (this->componentNorms == 0) {
//        this->allocComponentNorms();
//    }
//    for (int i = 0; i < this->getTDim(); i++) {
//        this->calcComponentNorm(i);
//    }
}

/** Calculate the norm of one component (NOT the squared norm!). */
template<int D>
void ProjectedNode<D>::calcComponentNorm(int i) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->componentNorms != 0);
//    VectorXd &c = this->getCoefs();
//    int kp1_d = this->getKp1_d();
//    this->componentNorms[i] = c.segment(i*kp1_d, kp1_d).norm();
}

/** Clear coefficients of generated nodes.
  *
  * The node structure is kept, only the coefficients are cleared. */
template<int D>
void ProjectedNode<D>::clearGenerated() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        for (int i = 0; i < this->getTDim(); i++) {
//            if (this->children[i] != 0) {
//                this->getFuncChild(i).clearGenerated();
//            }
//        }
//    }
}

/** Remove all generated nodes recursively
  *
  * All generated nodes are removed (nothing is kept) */
template<int D>
void ProjectedNode<D>::purgeGenerated() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        if (this->isEndNode()) {
//            this->deleteChildren();
//        } else {
//            for (int i = 0; i < this->getTDim(); i++) {
//                if (this->children[i] != 0) {
//                    this->getFuncChild(i).purgeGenerated();
//                }
//            }
//        }
//    }
}

/** Copies the coefs from the first ancestor that HAS coefs.
  *
  * Virtually resolved to stop recursion when a ProjectedNode is encountered.
  * In that case it is assumed that the genRootNode has been communicated in
  * MPI, otherwise you die. Returns the depth of the node it copied from. */
template<int D>
int ProjectedNode<D>::getGenRootCoefs(VectorXd &c) {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->hasCoefs()) {
//        MSG_FATAL("If you're running MPI, you need to communicate the " <<
//                  "genRootNodes prior to this function call. If you are " <<
//                  "NOT running MPI, something is wrong")
//    }
//    c = VectorXd::Zero(this->getNCoefs());
//    c = this->getCoefs();

//    return this->getDepth();
}

template class ProjectedNode<1> ;
template class ProjectedNode<2> ;
template class ProjectedNode<3> ;

