#ifndef FUNCTIONNODE_H
#define FUNCTIONNODE_H

#include "MWNode.h"

template<int D>
class FunctionNode : public MWNode<D> {
public:
    FunctionNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    FunctionNode(FunctionNode<D> &p, int cIdx);
    FunctionNode(const FunctionNode<D> &n);
    FunctionNode& operator=(const FunctionNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~FunctionNode() { }

    virtual double evalf(const double *r);
    virtual void clearGenerated() { NOT_IMPLEMENTED_ABORT; }
    virtual void purgeGenerated() { NOT_IMPLEMENTED_ABORT; }

    double integrate();
    double dotScaling(FunctionNode<D> &inpNode);
    double dotWavelet(FunctionNode<D> &inpNode);

    FunctionTree<D> &getFuncTree() { return static_cast<FunctionTree<D> &>(*this->tree); }
    FunctionNode<D> &getFuncParent() { return static_cast<FunctionNode<D> &>(*this->parent); }
    FunctionNode<D> &getFuncChild(int i) { return static_cast<FunctionNode<D> &>(this->children.getNode(i)); }

    const FunctionTree<D> &getFuncTree() const { return static_cast<const FunctionTree<D> &>(*this->tree); }
    const FunctionNode<D> &getFuncParent() const { return static_cast<const FunctionNode<D> &>(*this->parent); }
    const FunctionNode<D> &getFuncChild(int i) const { return static_cast<const FunctionNode<D> &>(this->children.getNode(i)); }

protected:
    double integrateLegendre();
    double integrateInterpolating();

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<MWNode<D> >(*this);
    }
};

#endif // FUNCTIONNODE_H
