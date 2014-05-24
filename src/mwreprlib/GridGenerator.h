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

#include "mwrepr_declarations.h"

template<int D>
class GridGenerator {
public:
    GridGenerator(int u = 0) : uniform(u), grid (0) { }
    virtual ~GridGenerator();

    void setUniform(int u) { this->uniform = u; }
    void generateGrid(MRGrid<D> &outGrid);

protected:
    int uniform;
    MRGrid<D> *grid;

    void splitNodeTable(GridNodeVector &nodeTable);
    virtual bool splitCheck(GridNode<D> *node);
};

#endif /* GRID_GENERATOR_H_ */
