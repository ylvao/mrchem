/**
*
*
*  \date Aug 14, 2009
*  \author Jonas Juselius <jonas.juselius@uit.no> \n
*  CTCC, University of Tromsø
*
*  Basic class for representing functions in a multiwavelet
* representation.
*/

#ifndef FUNCTIONTREE_S_H_
#define FUNCTIONTREE_S_H_

template<int D> class MultiResolutionAnalysis;
template<int D> class ProjectedNode;
template<int D> class FunctionTree;
template<int D> class FunctionNode;

template<int D>
class FunctionTree_S  {
public:
    FunctionTree_S(const MultiResolutionAnalysis<D> &mra, int maxNumberOfNodes);
    FunctionTree_S(const FunctionTree_S<D> &tree);
    FunctionTree_S<D> &operator=(const FunctionTree_S<D> &tree);
    virtual ~FunctionTree_S();

    void clear();

    double integrate() const;
    virtual double dot(const FunctionTree_S<D> &ket);
    virtual double evalf(const double *r);

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void orthogonalize(const FunctionTree_S<D> &tree);

    FunctionTree_S<D>& operator *=(double c);
    FunctionTree_S<D>& operator *=(const FunctionTree_S<D> &tree);
    FunctionTree_S<D>& operator +=(const FunctionTree_S<D> &tree);
    FunctionTree_S<D>& operator -=(const FunctionTree_S<D> &tree);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->funcTree_p->getEndNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->funcTree_p->rootBox.getNode(i)); }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->funcTree_p->getEndNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->funcTree_p->rootBox.getNode(i)); }

protected:
    // Static default parameters
    const static int tDim = (1 << D);

    int sizeTreeMeta; //The first part of the Tree is filled with metadata; reserved size:
    int sizeNode;     //The dynamical part of the tree is filled with nodes of size:
    int maxNodes;     //max number of nodes that can be defined
    int nNodes;       //number of nodes already defined

    double* tree_S_array; //Tree is defined as array of doubles, because C++ does not like void malloc

    FunctionTree<D>* funcTree_p;
    ProjectedNode<D>* lastNode;//pointer to the last active node

    ProjectedNode<D>* allocNodes(int Nalloc);
};

#endif /* FUNCTIONTREE_S_H_*/
