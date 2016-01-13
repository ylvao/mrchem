/**
 *  Simple n-dimensional node
 *
 *  Created on: May 29, 2009
 *      Author: jonas
 */

#ifndef MWNODE_H_
#define MWNODE_H_

#include "MRNode.h"
#include "ProjectedNode.h"

template<int D>
class MWNode : public MRNode<D> {
public:
    MWNode(MWTree<D> &t, const NodeIndex<D> &nIdx);
    MWNode(MWNode<D> &p, int cIdx);
    MWNode(const MWNode<D> &n);
    MWNode& operator=(const MWNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~MWNode();

    double estimateError(bool absPrec);

    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }

    double getSquareNorm() const { return this->squareNorm; }
    double getScalingNorm() const { return calcScalingNorm(); }
    double getWaveletNorm() const { return calcWaveletNorm(); }
    double getComponentNorm(int i) const { return this->componentNorms[i]; }
    bool hasComponentNorms() const;

    int getNCoefs() const { return this->coefs->size(); }
    virtual Eigen::VectorXd &getCoefs() { if (not this->isAllocated()) allocCoefs(); return *this->coefs; }
    virtual const Eigen::VectorXd &getCoefs() const { return *this->coefs; }

    virtual void setCoefs(const Eigen::VectorXd &c);
    virtual void zeroCoefs();

    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    MWTree<D>& getMWTree() { return static_cast<MWTree<D> &>(*this->tree); }
    MWNode<D>& getMWParent() { return static_cast<MWNode<D> &>(*this->parent); }
    MWNode<D>& getMWChild(int i) { return static_cast<MWNode<D> &>(*this->children[i]); }

    const MWTree<D>& getMWTree() const { return static_cast<const MWTree<D> &>(*this->tree); }
    const MWNode<D>& getMWParent() const { return static_cast<const MWNode<D> &>(*this->parent); }
    const MWNode<D>& getMWChild(int i) const { return static_cast<const MWNode<D> &>(*this->children[i]); }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const MWNode<T> &nd);

protected:
    double squareNorm;
    double *componentNorms; ///< 2^D components
    Eigen::VectorXd *coefs;

    virtual void allocCoefs(int nCoefs = -1);
    virtual void freeCoefs();

    void calcNorms();
    void zeroNorms();
    void clearNorms();

    double calcSquareNorm() const;
    double calcScalingNorm() const;
    virtual double calcWaveletNorm() const;

    void calcComponentNorms();
    virtual double calcComponentNorm(int i) const;

    inline void allocComponentNorms();
    inline void freeComponentNorms();

    void giveChildrenCoefs(bool overwrite = true);
    void copyCoefsFromChildren(Eigen::VectorXd &c);

    bool crop(double prec, NodeIndexSet *cropIdx = 0);
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
void MWNode<D>::freeComponentNorms() {
    if (this->componentNorms != 0) {
        delete [] this->componentNorms;
        this->componentNorms = 0;
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
