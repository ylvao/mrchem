/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#ifndef MRGRID_H_
#define MRGRID_H_

#include <Eigen/Core>
#include "mwrepr_declarations.h"

template<int D>
class MRGrid {
public:
    MRGrid(int order = -1, const GridNodeBox<D> *box = 0);
    virtual ~MRGrid();

    int getOrder() const { return this->order; }
    int getKp1() const { return this->kp1; }
    int getKp1_d() const { return this->kp1_d; }
    int getMaxScale() const { return this->maxScale; }
    int getMaxDepth() const { return this->maxDepth; }
    int getRootScale() const { return this->rootBox->getRootScale(); }

    void incrementNodeCount(int scale);
    void decrementNodeCount(int scale);
    int getNodesAtDepth(int depth) const;

    void clearEndNodeTable();
    void copyEndNodeTable(GridNodeVector &outTable);
    void resetEndNodeTable();

    GridNodeBox<D> &getRootBox() { return *this->rootBox; }
    const GridNodeBox<D> &getRootBox() const { return *this->rootBox; }

    const double *getLowerBounds() const { return this->rootBox->getLowerBounds(); }
    const double *getUpperBounds() const { return this->rootBox->getUpperBounds(); }

    int getNNodes(int depth = -1) const;
    int countBranchNodes(int depth = -1);
    int countLeafNodes(int depth = -1);
    int countQuadPoints(int depth = -1);

    void getQuadraturePoints(Eigen::MatrixXd &roots, const NodeIndex<D> &idx) const;
    void getQuadratureWeights(Eigen::VectorXd &weights, const NodeIndex<D> &idx) const;

    void getQuadraturePoints(Eigen::MatrixXd &roots) const;
    void getQuadratureWeights(Eigen::VectorXd &weights) const;

    void saveGrid(const std::string &file);
    void loadGrid(const std::string &file);

    friend class GridGenerator<D>;

protected:
    const static int tDim = (1 << D);
    int order;		
    int kp1;
    int kp1_d;
    int maxScale;
    int maxDepth;
    GridNodeBox<D> *rootBox;

    GridNodeVector endNodeTable;    ///< List of leaf nodes
    std::vector<int> nodesAtDepth;  ///< Number of nodes per scale

    void initializeRootNodes();
    void yieldChildren(GridNodeVector &nodeTable, const NodeIndexSet &idxSet);
};

#endif /* MRGRID_H_*/
