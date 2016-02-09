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

#ifndef FUNCTIONTREE_H_
#define FUNCTIONTREE_H_

#include "ProjectedNode.h"

template<int D>
class FunctionTree_S  {
 public:
    FunctionTree_S(const MultiResolutionAnalysis<D> &mra, int MaxNumberofNodes);
    FunctionTree_S(const FunctionTree_S<D> &tree);
    FunctionTree_S<D> &operator=(const FunctionTree_S<D> &tree);
    virtual ~FunctionTree_S();

  //The first part of the Tree is filled with metadata; reserved size:
  int SizeTreeMeta;
  //The dynamical part of the tree is filled with nodes of size:
  int SizeNode;

  //Tree is defined as array of doubles, because C++ does not like void malloc
  double* Tree_S_array;

  FunctionTree<D>* FunctionTree_p;

  int MaxNnodes;//max number of nodes that can be defined
  int Nnodes;//number of nodes already defined
  
  ProjectedNode<D>* LastNode;//pointer to the last active node

  ProjectedNode<D>* AllocNodes(int Nalloc);
  
  // Static default parameters
  const static int tDim = (1 << D);

    void clear();

    double integrate() const;
    virtual double dot(const FunctionTree_S<D> &ket);
    virtual double evalf(const double *r);

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void orthogonalize(const FunctionTree_S<D> &tree);
    void map(const RepresentableFunction<1> &func);

    FunctionTree_S<D>& operator *=(double c);
    FunctionTree_S<D>& operator *=(const FunctionTree_S<D> &tree);
    FunctionTree_S<D>& operator +=(const FunctionTree_S<D> &tree);
    FunctionTree_S<D>& operator -=(const FunctionTree_S<D> &tree);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->FunctionTree_p->getEndNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->FunctionTree_p->rootBox.getNode(i)); }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->FunctionTree_p->getEndNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->FunctionTree_p->rootBox.getNode(i)); }

    template<int T>
    friend std::ostream& operator <<(std::ostream &o, FunctionTree_S<T> &tree);

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT;
    }
};

template<int D>
std::ostream& operator<<(std::ostream &o, FunctionTree_S<D> &tree) {
    o << std::endl << "*FunctionTree_S: " << tree.name << std::endl;
    o << "  square norm: " << tree.squareNorm << std::endl;
    o << "  root scale: " << tree.getRootScale() << std::endl;
    o << "  order: " << tree.order << std::endl;
    o << "  nodes: " << tree.getNNodes() << std::endl;
    o << "  local end nodes: " << tree.endNodeTable.size() << std::endl;
    o << "  genNodes: " << tree.getNGenNodes() <<
            " (" << tree.getNAllocGenNodes() << ")" <<std::endl;
    o << "  nodes per scale: " << std::endl;
    for (int i = 0; i < tree.nodesAtDepth.size(); i++) {
        o << "    scale=" << i + tree.getRootScale() << "  nodes="
          << tree.nodesAtDepth[i] << std::endl;
    }
    return o;
}


#endif /* FUNCTIONTREE_H_*/
