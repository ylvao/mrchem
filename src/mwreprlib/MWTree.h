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

class Filter;
class ScalingBasis;

template<int D>
class MWTree : public MRTree<D> {
public:
    MWTree(const BoundingBox<D> *box, int k, int type);
    MWTree(const MRGrid<D> &grid, int type);
    MWTree(const MWTree<D> &tree);
    MWTree<D> &operator=(const MWTree<D> &tree);
    virtual ~MWTree();

    void setZero(bool clearTreeNorm = true);
    void setAutoClean(bool clean = true) { this->autoCleanGenerated = clean; }

    int getScalingType() const { return this->scalingType; }
    bool getAutoClean() const { return this->autoCleanGenerated; }
    double getSquareNorm() const { return this->squareNorm; }
    double estimateError(bool absPrec);

    const Filter &getFilter() { return *this->filter; }
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

    static void setDefaultSplitType(int type);
    static void setDefaultScalingType(int type);

protected:
    int scalingType;
    bool autoCleanGenerated;
    double squareNorm;

    Filter *filter;
    ScalingBasis *scalingFunc;

    Eigen::MatrixXd **tmpCoefs;   ///< temp memory
    Eigen::VectorXd **tmpVector;  ///< temp memory
    Eigen::VectorXd **tmpMWCoefs; ///< temp memory

    static int defaultSplitType;
    static int defaultScalingType;

    void setupFilters(int type);
    void setupScalingBasis(int type);

    void allocWorkMemory();
    void freeWorkMemory();

    void calcTreeNorm(MRNodeVector *work = 0);

//    void updateMissingScalingPart(const MWNodeVector &nodeList);
//    void setupWorkTable(MWNodeVector &wt);
//    void updateWorkTable(MWNodeVector &workTable);
//    void calcWorkTable(RepresentableObject<D> &func, MWNodeVector &workTable);
//    MWNodeVector &splitWorkNodes(RepresentableObject<D> &func, MWNodeVector &workTable);
//    MWNodeVector &seedWorkTable(RepresentableObject<D> &func, MWNodeVector &workTable, bool filter = false);

//    void seedTree(RepresentableObject<D> &func, bool _filter = false);
//    void calcTree(RepresentableObject<D> &func);
//    void growTree(RepresentableObject<D> &func, int maxIter = -1);

//    template<typename Tree>
//    std::vector<Tree> &copyVector(std::vector<Tree> &dest,
//                                  const std::vector<Tree> &src) {
//        dest.clear();
//        dest.insert(dest.begin(), src.begin(), src.end());
//        return dest;
//    }

//    template<typename Tree>
//    std::vector<Tree> &appendVector(std::vector<Tree> &dest,
//                                    const std::vector<Tree> &src) {
//        dest.insert(dest.end(), src.begin(), src.end());
//        return dest;
//    }

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
