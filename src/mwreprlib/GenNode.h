
/*
 *
 *  \date Oct 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef GENNODE_H_
#define GENNODE_H_

#include "FunctionNode.h"

template<int D> class ProjectedNode;

template<int D>
class GenNode: public FunctionNode<D> {
public:
    GenNode(FunctionNode<D> &p, int cIdx);
    GenNode(GenNode<D> &n);
    virtual ~GenNode();

    double evalf(const double *r);
    double getComponentNorm(int i);

    void clearGenerated();
    void purgeGenerated();

    void setCoefs(const Eigen::VectorXd &c);
    Eigen::VectorXd& getCoefs();
    const Eigen::VectorXd& getCoefs() const;

    void mwTransform(int kind);
    void cvTransform(int kind);

    const ProjectedNode<D> *getGenRootNode() const { return this->genRootNode; }
    ProjectedNode<D> *getGenRootNode() { return this->genRootNode; }

protected:
    void calcComponentNorms() { }
    inline double calcSquareNorm();
    inline double calcScalingNorm();
    inline double calcWaveletNorm() { return 0.0; }

    MWNode<D> *retrieveNode(int n, const double *r);
    MWNode<D> *retrieveNode(const NodeIndex<D> &idx, bool empty = false);

private:
    ProjectedNode<D> *genRootNode;

    void createChild(int i);
    void genChild(int i);

    void lockSiblings();
    void unlockSiblings();

    void regenerateCoefs();
    void releaseCoefs();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};


#endif /* GENNODE_H_ */
