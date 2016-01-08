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
    NOT_IMPLEMENTED_ABORT;
//    this->setIsEndNode();
//    if (not this->isForeign()) {
//        this->allocCoefs();
//        this->zeroCoefs();
//        this->zeroNorms();
//    }
}

template<int D>
ProjectedNode<D>::ProjectedNode(const ProjectedNode<D> &n)
        : FunctionNode<D>(n) {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->isForeign()) {
//        this->allocCoefs();
//        this->zeroCoefs();
//        this->zeroNorms();
//    }
}

/* Recurcive node constructor*/
template<int D>
void ProjectedNode<D>::copyChildren(const MRNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
//    if (node.isBranchNode()) {
//        this->allocKindergarten();
//        this->setIsBranchNode();
//        this->clearIsEndNode();
//    }
//    int myRank = this->getRankId();
//    for (int cIdx = 0; cIdx < node.getNChildren(); cIdx++) {
//        const MRNode<D> &yourChild = node.getMRChild(cIdx);
//        int childRank = yourChild.getRankId();
//        this->setRankId(childRank); //Rank is copied from parent
//        ProjectedNode<D> *myChild = new ProjectedNode(*this, cIdx);
//        this->setRankId(myRank);
//        myChild->copyChildren(yourChild);
//        this->children[cIdx] = myChild;
//    }
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty ProjectedNode child node
  * with the appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChild(int cIdx) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->children != 0);
//    assert(this->children[cIdx] == 0);
//    ProjectedNode<D> *child = new ProjectedNode<D>(*this, cIdx);
//    this->children[cIdx] = child;
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
double ProjectedNode<D>::calcComponentNorm(int i) {
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

template class ProjectedNode<1> ;
template class ProjectedNode<2> ;
template class ProjectedNode<3> ;

