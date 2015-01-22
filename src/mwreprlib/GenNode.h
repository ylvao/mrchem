
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
    GenNode();
    GenNode(FunctionNode<D> *p, int cIdx);
    GenNode(const GenNode<D> &nd, FunctionTree<D> *t);
    GenNode(const GenNode<D> &nd, FunctionNode<D> *p);
    GenNode<D> &operator=(const GenNode<D> &nd);
    virtual ~GenNode();

    double evalf(const double *r);

    void clearCoefs();
    Eigen::VectorXd &getCoefs();
    const Eigen::VectorXd &getCoefs() const;
    void setCoefs(const Eigen::VectorXd &c);

    void mwTransform(int kind);
    void cvTransform(int kind);

    int getGenRootCoefs(Eigen::VectorXd &c);
    void clearGenerated();
    void purgeGenerated();

    const ProjectedNode<D> *getGenRootNode() const { return this->genRootNode; }
    ProjectedNode<D> *getGenRootNode() { return this->genRootNode; }

    void incrementNodeWeight(int i, double w);
    double getComponentNorm(int i);

    /** Calculate GenNode component norms, ie. do nothing.
    We only have one component norm == squareNorm. The component norm is
    calculated lazily by getComponentNorm(), since we do not necessarily have
    any coefs yet. */
    void calcComponentNorms() { }
    bool splitCheck(double prec = -1.0) { return false; }
    void allocCoefs(int nCoefs = -1);

protected:
    void createChildren();
    void createChild(int i);
    void genChildren(bool genEmpty = false);
    void giveChildrenScaling(bool overwrite = true);
    void giveSiblingsScaling(Eigen::VectorXd &coefs);

    /** Calculate GenNode component norms: ie. do nothing. */
    void calcComponentNorm(int i) {}
    inline double calcSquareNorm();
    inline double calcScalingNorm();
    inline double calcWaveletNorm() { return 0.0; }

    Eigen::VectorXd &getCoefsNoLock();
    void setCoefsNoLock(const Eigen::VectorXd &c);

    MWNode<D> *retrieveNode(int n, const double *r);
    MWNode<D> *retrieveNode(const NodeIndex<D> &idx, bool genEmpty = false);

private:
    ProjectedNode<D> *genRootNode;

    void lockSiblings();
    void unlockSiblings();
    void regenerateCoefs();
    void releaseCoefs();
    void makeIndexList(int *idxList, int depth, const NodeIndex<D> &idx);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};


#endif /* GENNODE_H_ */
