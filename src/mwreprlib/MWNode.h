/**
 *  Simple n-dimensional node
 *
 *  Created on: May 29, 2009
 *      Author: jonas
 */

#ifndef MWNODE_H_
#define MWNODE_H_

#include <Eigen/Core>

#include "MRNode.h"

template<int D> class MWTree;

template<int D>
class MWNode : public MRNode<D> {
public:
    MWNode();
    MWNode(MRTree<D> &t, const NodeIndex<D> &nIdx);
    MWNode(MWNode<D> *p, int cIdx);
    MWNode(const MWNode<D> &nd, MWNode<D> *p);
    MWNode(const MWNode<D> &nd, MWTree<D> *t);
    MWNode<D> &operator=(const MWNode<D> &nd);
    virtual ~MWNode();

    inline bool hasComponentNorms() const;
    virtual double getComponentNorm(int i);
    inline double getSquareNorm() const;
    inline double getScalingNorm();
    inline double getWaveletNorm();

    void calcNorms();
    void clearNorms();
    void zeroNorms();

    virtual Eigen::VectorXd &getCoefs();
    virtual const Eigen::VectorXd &getCoefs() const;

    int getNCoefs() const { return this->coefs->size(); }
    virtual void setCoefs(const Eigen::VectorXd &c);
    virtual void freeCoefs();
    virtual void zeroCoefs();

    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    virtual bool splitCheck(double prec = -1.0) = 0;
    virtual void allocCoefs(int nCoefs = -1);
    virtual void calcComponentNorms() = 0;
    virtual void purgeGenerated() = 0;

    double estimateError(bool absPrec);
    double getNodeWeight(int i) { return this->weight[i]; }
    void clearNodeWeight(int i) { this->weight[i] = 0.0;}
    virtual void incrementNodeWeight(int i, double w = 1.0) {
        this->lockNode();
        this->weight[i] += w;
        this->unlockNode();
    }

    inline MWTree<D>& getMWTree() {
        return static_cast<MWTree<D> &>(*this->tree);
    }
    inline MWNode<D>& getMWParent() {
        return static_cast<MWNode<D> &>(*this->parent);
    }
    inline MWNode<D>& getMWChild(int i) {
        return static_cast<MWNode<D> &>(*this->children[i]);
    }
    inline const MWTree<D>& getMWTree() const {
        return static_cast<const MWTree<D> &>(*this->tree);
    }
    inline const MWNode<D>& getMWParent() const {
        return static_cast<const MWNode<D> &>(*this->parent);
    }
    inline const MWNode<D>& getMWChild(int i) const {
        return static_cast<const MWNode<D> &>(*this->children[i]);
    }

    friend class MWTree<D>;

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const MWNode<T> &nd);

protected:
    double weight[2];
    double squareNorm;
    double *componentNorms; ///< 2^D components

    Eigen::VectorXd *coefs;

    bool diffBranch(const MWNode<D> &rhs) const;
    bool diffCoefs(const MWNode<D> &rhs) const;
    bool diffNorms(const MWNode<D> &rhs) const;

    virtual double calcSquareNorm();
    virtual double calcScalingNorm();
    virtual double calcWaveletNorm();

    void getChildrenQuadRoots(std::vector<Eigen::MatrixXd *> &quadPts);
//    bool seedCheck(RepresentableObject<D> &func);

    void reallocCoefs(int nCoefs = -1);
    inline void allocComponentNorms();
    void freeComponentNorms();

    virtual void genChildren(bool genEmpty = false) = 0;
    virtual void giveChildrenScaling(bool overwrite = true) {}
    void copyScalingCoefsFromChildren(Eigen::VectorXd &scoefs);
    bool cropChildren(double prec, std::set<const NodeIndex<D> *,
                      NodeIndexComp<D> > *cropIdx = 0);

    void reCompress(bool overwrite = true);

    mpi::request isendCoefs(int who, int tag, int comp = -1);
    mpi::request ireceiveCoefs(int who, int tag, int comp = -1);

private:
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
        NOT_IMPLEMENTED_ABORT;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

template<int D>
double MWNode<D>::getComponentNorm(int i) {
    assert(i >= 0 and i < this->getTDim());
    assert(this->componentNorms != 0);
    assert(this->componentNorms[i] != -1.0);
    return this->componentNorms[i];
}

template<int D>
double MWNode<D>::getSquareNorm() const {
    assert(this->squareNorm >= 0.0);
    return this->squareNorm;
}

template<int D>
double MWNode<D>::getScalingNorm() {
    assert(this->componentNorms != 0);
    return getComponentNorm(0);
}

template<int D>
double MWNode<D>::getWaveletNorm() {
    return calcWaveletNorm();
}

template<int D>
bool MWNode<D>::hasComponentNorms() const {
    if (this->componentNorms != 0) {
        if (this->componentNorms[0] == -1.0) {
            return false;
        }
        return true;
    }
    return false;
}

template<int D>
void MWNode<D>::allocComponentNorms() {
    if (this->componentNorms == 0) {
        this->componentNorms = new double[this->getTDim()];
        for (int i = 0; i < this->getTDim(); i++) {
            this->componentNorms[i] = -1.0;
        }
    }
}

template<int D>
std::ostream& operator<<(std::ostream &o, const MWNode<D> &nd) {
    MRNode<D>::operator<<(o, nd);
    o << " sqNorm=" << nd.squareNorm;
    if (nd.hasCoefs()) {
        o << " Coefs={";
        o << nd.getCoefs()[0] << ", " <<
                nd.getCoefs()[nd.getNCoefs() - 1] << "}";
    }
    return o;
}

#endif /* MWNODE_H_ */
