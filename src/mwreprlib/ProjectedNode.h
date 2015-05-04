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
    ProjectedNode();
    ProjectedNode(FunctionTree<D> &t, const GridNode<D> &gNode);
    ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    ProjectedNode(ProjectedNode<D> &p, int cIdx);
    ProjectedNode(ProjectedNode<D> &p, int cIdx, const GridNode<D> &gNode);
    ProjectedNode(const ProjectedNode<D> &nd, ProjectedNode<D> *p);
    ProjectedNode(const ProjectedNode<D> &nd, FunctionTree<D> *t);
    ProjectedNode<D> &operator=(const ProjectedNode<D> &nd);
    virtual ~ProjectedNode() { }

    int getGenRootCoefs(Eigen::VectorXd &c);
    void clearGenerated();
    void purgeGenerated();

    bool splitCheck(double prec = -1.0);
    void calcComponentNorms();

    void genChildren(bool genEmpty = false);
    void createChildren();
    void deleteChildren() {
        MWNode<D>::deleteChildren();
        this->setIsEndNode();
    }

protected:
    void copyRecursive(GridNode<D> &gNode);
    void calcComponentNorm(int i);

//    void copyScalingCoefs(const ProjectedNode<D> &nd);
//    void giveChildrenScaling(bool overwrite = true);

private:
    void createChild(int i);
    void genChild(int i);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};

#endif /* PROJECTEDNODE_H_ */
