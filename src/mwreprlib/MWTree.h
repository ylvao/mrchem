/**
 *
 * \date Jun 5, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 */

#ifndef MWTREE_H_
#define MWTREE_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>
#include <Eigen/Core>

#include "MRTree.h"
#include "MWNode.h"
#include "MultiResolutionAnalysis.h"

template<int D>
class MWTree : public MRTree<D> {
public:
    MWTree(const MultiResolutionAnalysis<D> &mra);
    MWTree(const MWTree<D> &tree);
    virtual ~MWTree();

    void setZero();

    double estimateError(bool absPrec);
    double getSquareNorm() const { return this->squareNorm; }

    int getOrder() const { return this->order; }
    int getKp1() const { return this->order + 1; }
    int getKp1_d() const { return this->kp1_d; }
    const MultiResolutionAnalysis<D> &getMRA() const { return this->MRA; }

    void crop(double thrs = -1.0, bool absPrec = true);
    void mwTransform(int type, bool overwrite = true);

    MWNode<D> &getEndMWNode(int i) { return static_cast<MWNode<D> &>(this->getEndNode(i)); }
    MWNode<D> &getRootMWNode(int i) { return static_cast<MWNode<D> &>(this->rootBox.getNode(i)); }

    const MWNode<D> &getEndMWNode(int i) const { return static_cast<const MWNode<D> &>(this->getEndNode(i)); }
    const MWNode<D> &getRootMWNode(int i) const { return static_cast<const MWNode<D> &>(this->rootBox.getNode(i)); }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, MWTree<T> &tree);
    friend class MWNode<D>;

protected:
    const MultiResolutionAnalysis<D> MRA;

    // Constant parameters that are derived internally
    const int order;
    const int kp1_d;

    // Tree data
    double squareNorm;

    Eigen::MatrixXd **tmpCoefs;   ///< temp memory
    Eigen::VectorXd **tmpVector;  ///< temp memory
    Eigen::VectorXd **tmpMWCoefs; ///< temp memory

    void allocWorkMemory();
    void freeWorkMemory();

    inline Eigen::MatrixXd &getTmpScalingCoefs();
    inline Eigen::VectorXd &getTmpScalingVector();
    inline Eigen::VectorXd &getTmpMWCoefs();

    void calcSquareNorm(const MRNodeVector *work = 0);

    void mwTransformDown(bool overwrite = true);
    void mwTransformUp(bool overwrite = true);

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
        NOT_IMPLEMENTED_ABORT;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT;
//        setupFilters(scalingType);
//        setupScalingBasis(scalingType);
//        freeWorkMemory();
//        allocWorkMemory();

//        if (this->isScattered()) {
//            if (node_group.size() < 2) {
//                MSG_WARN("Reading distributed tree in serial. " <<
//                         "Tree is incomplete and unpure.");
//                resetEndNodeTable();
//                return;
//            }
//        }
//        mwTransformUp();
//        resetEndNodeTable();
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

template<int D>
Eigen::MatrixXd& MWTree<D>::getTmpScalingCoefs() {
    int thread = omp_get_thread_num();
    return *this->tmpCoefs[thread];
}

template<int D>
Eigen::VectorXd& MWTree<D>::getTmpScalingVector() {
    int thread = omp_get_thread_num();
    return *this->tmpVector[thread];
}

template<int D>
Eigen::VectorXd& MWTree<D>::getTmpMWCoefs() {
    int thread = omp_get_thread_num();
    return *this->tmpMWCoefs[thread];
}

template<int T>
std::ostream& operator<<(std::ostream &o, MWTree<T> &tree) {
    o << "*MWTree: " << tree.name << std::endl;
    o << "  order: " << tree.order << std::endl;
    o << "  nodes: " << tree.nNodes << std::endl;
    o << "  genNodes: " << tree.getNGenNodes() <<
            " (" << tree.getNAllocGenNodes() << ")" << std::endl;
    o << "  nodes per scale: " << std::endl;
    return o;
}

#endif /* MWTREE_H_ */
