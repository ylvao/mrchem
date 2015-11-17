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

class MWFilter;
class ScalingBasis;
template<int D> class MWProjector;

template<int D>
class MWTree : public MRTree<D> {
public:
    MWTree(const BoundingBox<D> &box, int k, int type);
    MWTree(const MRGrid<D> &tree, int type);
    MWTree(const MWTree<D> &tree);
    virtual ~MWTree();

    void setZero();
    virtual void clear() = 0;

    int getScalingType() const { return this->scalingType; }
    double getSquareNorm() const { return this->squareNorm; }
    double estimateError(bool absPrec);

    const MWFilter &getFilter() { return *this->filter; }
    const ScalingBasis &getScalingFunctions() const { return *this->scalingFunc; }

    inline Eigen::MatrixXd &getTmpScalingCoefs();
    inline Eigen::VectorXd &getTmpScalingVector();
    inline Eigen::VectorXd &getTmpMWCoefs();

    void mwTransformDown(bool overwrite = true);
    void mwTransformUp(bool overwrite = true);

    void crop(double thrs = -1.0, bool absPrec = true);

    MWNode<D>& getRootMWNode(int rIdx);
    MWNode<D>& getRootMWNode(const NodeIndex<D> &nIdx);

    const MWNode<D>& getRootMWNode(int rIdx) const;
    const MWNode<D>& getRootMWNode(const NodeIndex<D> &nIdx) const;

    MWNode<D>& getEndMWNode(int i);
    const MWNode<D>& getEndMWNode(int i) const;

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, MWTree<T> &tree);
    friend class MWProjector<D>;

protected:
    int scalingType;
    double squareNorm;

    MWFilter *filter;
    ScalingBasis *scalingFunc;

    Eigen::MatrixXd **tmpCoefs;   ///< temp memory
    Eigen::VectorXd **tmpVector;  ///< temp memory
    Eigen::VectorXd **tmpMWCoefs; ///< temp memory

    void setupFilters(int type, int k);
    void setupScalingBasis(int type, int k);

    void allocWorkMemory();
    void freeWorkMemory();

    double calcSquareNorm(MRNodeVector *work = 0);

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
    o << "  scaling type: " << tree.scalingType << std::endl;
    o << "  order: " << tree.order << std::endl;
    o << "  nodes: " << tree.nNodes << std::endl;
    o << "  genNodes: " << tree.getNGenNodes() <<
            " (" << tree.getNAllocGenNodes() << ")" << std::endl;
    o << "  nodes per scale: " << std::endl;
    return o;
}

#endif /* MWTREE_H_ */
