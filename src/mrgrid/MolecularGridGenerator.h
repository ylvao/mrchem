/*
 *
 *  \date May 23, 2014
 *  \author Stig Rune Jensen
 *          CTCC, University of Tromsø
 *
 *  Base class to generate initial grids
 */

#ifndef MOLECULAR_GRID_GENERATOR_H_
#define MOLECULAR_GRID_GENERATOR_H_

#include "GridGenerator.h"

class Molecule;

class MolecularGridGenerator : public GridGenerator<3> {
public:
    MolecularGridGenerator() : depth(-1), width(-1), molecule(0) { }
    virtual ~MolecularGridGenerator();

    void setDepth(int d) { this->width = d; }
    void setWidth(int w) { this->width = w; }

    void generateGrid(MRGrid<3> &grid, Molecule &mol);

protected:
    int depth;
    int width;
    Molecule *molecule;

    bool splitCheck(GridNode<3> *node);
};

#endif /* GRID_GENERATOR_H_ */
