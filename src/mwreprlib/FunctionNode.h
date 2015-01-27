#ifndef FUNCTIONNODE_H
#define FUNCTIONNODE_H

#include "MWNode.h"

template<int D> class FunctionTree;

template<int D>
class FunctionNode : public MWNode<D> {
public:
    FunctionNode();
    FunctionNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    FunctionNode(FunctionNode<D> *p, int cIdx);
    FunctionNode(const FunctionNode<D> &nd, FunctionNode<D> *p);
    FunctionNode(const FunctionNode<D> &nd, FunctionTree<D> *t);
    FunctionNode<D> &operator=(const FunctionNode<D> &nd);
    virtual ~FunctionNode() { }

    virtual double evalf(const double *r);

    virtual int getGenRootCoefs(Eigen::VectorXd &c) = 0;
    virtual void clearGenerated() = 0;
    virtual void calcComponentNorms() = 0;

    void operator*=(double c);
    double integrate();
    double scalingInnerProduct(FunctionNode<D> &inpNode);
    double waveletInnerProduct(FunctionNode<D> &inpNode);

//    void addCoefs(Eigen::VectorXd &expCoefs, Eigen::MatrixXd &mwCoefs);
//    void addCoefs(double a, FunctionNode<D> &lhNode,
//                  double b, FunctionNode<D> &rhNode);

//    void multiplyCoefs(Eigen::VectorXd &expCoefs, Eigen::MatrixXd &mwCoefs);
//    void multiplyCoefs(double a, FunctionNode<D> &lhNode,
//                       double b, FunctionNode<D> &rhNode);

    inline FunctionTree<D> &getFuncTree();
    inline FunctionNode<D> &getFuncParent();
    inline FunctionNode<D> &getFuncChild(int cIdx);

    inline const FunctionTree<D> &getFuncTree() const;
    inline const FunctionNode<D> &getFuncParent() const;
    inline const FunctionNode<D> &getFuncChild(int cIdx) const;

    friend class FunctionTree<D>;

protected:
    double integrateLegendre();
    double integrateInterpolating();

    MRNode<D> *retrieveNode(int n, const double *r);
    MRNode<D> *retrieveNode(const NodeIndex<D> &idx);

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<MWNode<D> >(*this);
    }
};

/** Static cast of MWTree reference to FunctionTree reference.
  *
  * FunctionTrees only contain FunctionNodes, and FunctionNodes are always part
  * of a FunctionTree. Still the tree pointer of the node is defined in the
  * MWNode base class, and thus its tree pointer must be of MWTree type. This
  * routine returns the tree as a FunctionTree. */
template<int D>
const FunctionTree<D> &FunctionNode<D>::getFuncTree() const {
    assert(this->tree != 0);
    return static_cast<const FunctionTree<D> &>(*this->tree);
}

template<int D>
FunctionTree<D> &FunctionNode<D>::getFuncTree() {
    assert(this->tree != 0);
    return static_cast<FunctionTree<D> &>(*this->tree);
}

template<int D>
const FunctionNode<D> &FunctionNode<D>::getFuncChild(int cIdx) const {
    assert(this->children != 0);
    assert(this->children[cIdx] != 0);
    return static_cast<const FunctionNode<D> &>(*this->children[cIdx]);
}

template<int D>
FunctionNode<D> &FunctionNode<D>::getFuncChild(int cIdx) {
    assert(this->children != 0);
    assert(this->children[cIdx] != 0);
    return static_cast<FunctionNode<D> &>(*this->children[cIdx]);
}

template<int D>
const FunctionNode<D> &FunctionNode<D>::getFuncParent() const {
    assert(this->parent != 0);
    return static_cast<const FunctionNode<D> &>(*this->parent);
}

template<int D>
FunctionNode<D> &FunctionNode<D>::getFuncParent() {
    assert(this->parent != 0);
    return static_cast<FunctionNode<D> &>(*this->parent);
}

#endif // FUNCTIONNODE_H
