#ifndef MWREPR_DECLARATIONS_H_
#define MWREPR_DECLARATIONS_H_

#include <vector>
#include <set>

template <int D> class BoundingBox;
template <int D> class NodeBox;
template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

template <int D> class RepresentableFunction;

template <int D> class MRTree;
template <int D> class MWTree;
template <int D> class FunctionTree;

template <int D> class MRNode;
template <int D> class MWNode;
template <int D> class FunctionNode;
template <int D> class ProjectedNode;
template <int D> class GenNode;

template <int D> class TreeBuilder;
template <int D> class TreeAdaptor;
template <int D> class MWAdaptor;
template <int D> class TreeProjector;
template <int D> class FunctionProjector;


template <int D> class TreeIterator;
template <int D> class LebesgueIterator;
template <int D> class HilbertIterator;
template <int D> class HilbertPath;

#define MRNodeVector std::vector<MRNode<D> *>
#define MWNodeVector std::vector<MWNode<D> *>
#define NodeIndexSet std::set<const NodeIndex<D> *, NodeIndexComp<D> >

#endif /* MWREPR_DECLARATIONS_H_*/
