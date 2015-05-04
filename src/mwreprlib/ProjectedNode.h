/**
 *
 *  \date Aug 14, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 */

#ifndef PROJECTEDNODE_H_
#define PROJECTEDNODE_H_

#include "FunctionNode.h"

template<int D>
class ProjectedNode: public FunctionNode<D> {
public:
    ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    ProjectedNode(ProjectedNode<D> &p, int cIdx);
    ProjectedNode<D> &operator=(const ProjectedNode<D> &nd);
    virtual ~ProjectedNode() { }

    int getGenRootCoefs(Eigen::VectorXd &c);
    void clearGenerated();
    void purgeGenerated();

    void calcComponentNorms();

    void genChildren(bool genEmpty = false);
    void createChildren();
    void deleteChildren() {
        MWNode<D>::deleteChildren();
        this->setIsEndNode();
    }

private:
    void calcComponentNorm(int i);
    void createChild(int i);
    void genChild(int i);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};

#endif /* PROJECTEDNODE_H_ */
