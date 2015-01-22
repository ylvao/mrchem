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
    ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    ProjectedNode(ProjectedNode<D> *p, int cIdx);
    ProjectedNode(const ProjectedNode<D> &nd, ProjectedNode<D> *p);
    ProjectedNode(const ProjectedNode<D> &nd, FunctionTree<D> *t);
    ProjectedNode<D> &operator=(const ProjectedNode<D> &nd);
    virtual ~ProjectedNode() { }

    int getGenRootCoefs(Eigen::VectorXd &c);
    void clearGenerated();
    void purgeGenerated();

    bool splitCheck(double prec = -1.0);
    void calcComponentNorms();

protected:
    void createChildren();
    void deleteChildren() {
        MWNode<D>::deleteChildren();
        this->setIsEndNode();
    }

    void genChildren(bool genEmpty = false);
    void genChild(int i);

    void calcComponentNorm(int i);
    double getWaveletThreshold(int factor, double prec);

//    void copyScalingCoefs(const ProjectedNode<D> &nd);
//    void giveChildrenScaling(bool overwrite = true);

private:
    void createChild(int i);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};

#endif /* PROJECTEDNODE_H_ */
