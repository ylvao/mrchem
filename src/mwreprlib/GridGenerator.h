/*
 *
 *  \date May 23, 2014
 *  \author Stig Rune Jensen
 *          CTCC, University of Tromsø
 *
 *  Base class to generate initial grids
 */

#ifndef GRID_GENERATOR_H_
#define GRID_GENERATOR_H_

#include "constants.h"
#include "mwrepr_declarations.h"

template<int D>
class GridGenerator {
public:
    GridGenerator(int uScale = MinScale);
    virtual ~GridGenerator();

    void setUniformScale(int u) { this->uniformScale = u; }
    void generateGrid(MRGrid<D> &outGrid);
protected:
    int uniformScale;

    void initGrid(MRGrid<D> &outGrid);
    void clearGrid();
    void buildGrid();

    virtual bool splitCheck(const GridNode<D> *node);

    const MRGrid<D> *getGrid() const { return this->grid; }
private:
    MRGrid<D> *grid;
    int quadOrder;
    int uniformDepth;

    void splitNodeTable(GridNodeVector &nodeTable);
    bool updateNodeTable(GridNodeVector &nodeTable, NodeIndexSet &idxSet);
};

#endif /* GRID_GENERATOR_H_ */
