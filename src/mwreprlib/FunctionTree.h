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

#include "MWTree.h"
#include "RepresentableFunction.h"

template<int D> class FunctionNode;
template<int D> class BoundingBox;
template<int D> class NodeIndex;
template<int D> class MWProjector;

template<int D>
class FunctionTree: public MWTree<D>, public RepresentableFunction<D> {
public:
    FunctionTree(const BoundingBox<D> *box = 0, int k = -1,  int type = Interpol);
    FunctionTree(const MRGrid<D> &grid, int type = Interpol);
    FunctionTree(const FunctionTree<D> &tree);
    FunctionTree<D> &operator=(const FunctionTree<D> &tree);
    virtual ~FunctionTree();

    void clear();
    void clearGenNodes();
    void purgeGenNodes();

    double integrate();
    virtual double innerProduct(FunctionTree<D> &tree);
    virtual double evalf(const double *r) const;

    void plotCleanup() { this->purgeGenNodes(); }
    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    FunctionNode<D>& getRootFuncNode(int rIdx);
    FunctionNode<D>& getRootFuncNode(const NodeIndex<D> &nIdx);

    const FunctionNode<D>& getRootFuncNode(int rIdx) const;
    const FunctionNode<D>& getRootFuncNode(const NodeIndex<D> &nIdx) const;

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void orthogonalize(const FunctionTree<D> &tree);
    void map(const RepresentableFunction<1> &func);
    FunctionTree<D>& operator *=(double c);
    FunctionTree<D>& operator *=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator +=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator -=(const FunctionTree<D> &tree);

    template<int T>
    friend std::ostream& operator <<(std::ostream &o, FunctionTree<T> &tree);

    friend class MWProjector<D>;

protected:
    void initializeRootNodes();
    void initializeNodesRecursive(const MRGrid<D> &grid);

private:
//    HilbertPathTable<D> hilbertPathTable;

    void findMissingInnerProd(FunctionTree<D> &tree,
                              std::set<FunctionNode<D> *> &missing);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT;
    }
};

template<int D>
std::ostream& operator<<(std::ostream &o, FunctionTree<D> &tree) {
    o << std::endl << "*FunctionTree: " << tree.name << std::endl;
    o << "  squareNorm: " << tree.squareNorm << std::endl;
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
