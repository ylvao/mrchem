#ifndef MWREPR_DECLARATIONS_H_
#define MWREPR_DECLARATIONS_H_

#include <vector>
#include <set>

template<int D> class RepresentableFunction;

template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

template <int D> class GridGenerator;
template <int D> class GridNodeBox;

template <int D> class MRTree;
template <int D> class MRNode;

template <int D> class MWTree;
template <int D> class MWNode;

template <int D> class MRGrid;
template <int D> class GridNode;

template<int D> class MWAdaptor;
template<int D> class GridAdaptor;
template<int D> class FunctionTree;

#define MRNodeVector std::vector<MRNode<D> *>
#define MWNodeVector std::vector<MWNode<D> *>
#define NodeIndexSet std::set<const NodeIndex<D> *, NodeIndexComp<D> >

#endif /* MWREPR_DECLARATIONS_H_*/
