/*
 *
 *
 *  \date Oct 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "ProjectedNode.h"
#include "GenNode.h"

using namespace std;
using namespace Eigen;

template<int D>
GenNode<D>::GenNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(&p) {
    this->setIsGenNode();
}

template<int D>
GenNode<D>::GenNode(GenNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(p.genRootNode) {
    this->setIsGenNode();
}

template<int D>
GenNode<D>::GenNode(const GenNode<D> &n)
        : FunctionNode<D> (n),
          genRootNode(n.genRootNode) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
GenNode<D>::~GenNode() {
}

template<int D>
void GenNode<D>::createChild(int i) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MRNode<D> *child = new GenNode<D>(*this, cIdx);
    this->children[cIdx] = child;
}

template<int D>
void GenNode<D>::regenerateCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::releaseCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::setCoefs(const VectorXd &c) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
VectorXd& GenNode<D>::getCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
const VectorXd& GenNode<D>::getCoefs() const {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::clearGenerated() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::purgeGenerated() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::mwTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::cvTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::lockSiblings() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::unlockSiblings() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
double GenNode<D>::getComponentNorm(int i) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
double GenNode<D>::calcSquareNorm() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
double GenNode<D>::calcScalingNorm() {
    NOT_IMPLEMENTED_ABORT;
}

template class GenNode<1> ;
template class GenNode<2> ;
template class GenNode<3> ;
