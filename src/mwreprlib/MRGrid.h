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
#include "TelePrompter.h"

template<int D>
class MRGrid {
public:
    MRGrid(int order = -1, const GridNodeBox<D> *box = 0);
    virtual ~MRGrid();

    int getTDim() const { return this->tDim; }
    int getOrder() const { return this->order; }
    int getKp1() const { return this->kp1; }
    int getKp1_d() const { return this->kp1_d; }
    int getMaxScale() const { return this->maxScale; }
    int getMaxDepth() const { return this->maxDepth; }
    int getRootScale() const { return this->rootBox->getRootScale(); }
    const double *getOrigin() const { return this->rootBox->getOrigin(); }

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

    void getQuadPoints(Eigen::MatrixXd &points);
    void getQuadWeights(Eigen::VectorXd &weights);

    void saveGrid(const std::string &file);
    void loadGrid(const std::string &file);

    friend std::ostream& operator<<(std::ostream &o, MRGrid<D> &grid) {
	o << std::endl << std::endl;
	o << "===============================================" << std::endl;
	o << "|                   MRGrid                    |" << std::endl;
	o << "|---------------------------------------------|" << std::endl;
	int root = grid.getRootScale();
	int nScales = grid.nodesAtDepth.size();
	o << "| Scale "  << "   AllNodes " << "  LeafNodes " << "   QuadPoints |" << std::endl;
	for (int depth = 0; depth < nScales; depth++) {
	    int scale = depth + root;
	    o << "|" << std::setw(5) << scale << std::setw(12) << grid.getNNodes(depth);
	    o << std::setw(12) << grid.countLeafNodes(depth);
	    o << std::setw(14) << grid.countQuadPoints(depth);
	    o << "  |" << std::endl;
	}
	o << "|---------------------------------------------|" << std::endl;
	o << "| Total " << std::setw(10) << grid.getNNodes();
	o << std::setw(12) << grid.countLeafNodes();
	o << std::setw(14) << grid.countQuadPoints() << "  |" << std::endl;
	o << "===============================================" << std::endl << std::endl;
	return o;
    }

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
