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

#include "Eigen/Core"
#include "MRTree.h"

template<int D>
class MRGrid : public MRTree<D> {
public:
    MRGrid(const MRGrid<D> &grid);
    MRGrid(const BoundingBox<D> &box, int k);
    virtual ~MRGrid();
    void clear();

    int countQuadPoints(int depth = -1);
    void getQuadPoints(Eigen::MatrixXd &points);
    void getQuadWeights(Eigen::VectorXd &weights);

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    GridNode<D>& getRootGridNode(int rIdx);
    GridNode<D>& getRootGridNode(const NodeIndex<D> &nIdx);

    const GridNode<D>& getRootGridNode(int rIdx) const;
    const GridNode<D>& getRootGridNode(const NodeIndex<D> &nIdx) const;

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, MRGrid<T> &grid);

    friend class GridGenerator<D>;
protected:
    void initializeRootNodes();
};

template<int D>
std::ostream& operator<<(std::ostream &o, MRGrid<D> &grid) {
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

#endif /* MRGRID_H_*/
