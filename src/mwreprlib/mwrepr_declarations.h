#ifndef MWREPR_DECLARATIONS_H_
#define MWREPR_DECLARATIONS_H_

#include <vector>
#include <set>

class Filter;
class ScalingBasis;

template <int D> class TreeBuilder;
template <int D> class MWTree;
template <int D> class MWNode;
template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

template <int D> class FunctionTree;
template <int D> class FunctionTreeProxy;
template <int D> class FunctionNode;
template <int D> class ProjectedNode;
template <int D> class GenNode;

template <int D> class MRGrid;
template <int D> class GridNode;
template <int D> class GridNodeBox;

template <int D> class TreeIterator;
template <int D> class HilbertPathTable;
template <int D> class HilbertIterator;
template <int D> class LebesgueIterator;
template <int D> class BoundingBox;
template <int D> class NodeBox;

#define GridNodeVector std::vector<GridNode<D> *>
#define MWNodeVector std::vector<MWNode<D> *>
#define NodeIndexSet std::set<const NodeIndex<D> *, NodeIndexComp<D> > 

#endif /* MWREPR_DECLARATIONS_H_*/
