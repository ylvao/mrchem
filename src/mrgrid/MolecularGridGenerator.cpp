/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include "MolecularGridGenerator.h"
#include "Molecule.h"
#include "MRGrid.h"


MolecularGridGenerator::~MolecularGridGenerator() {
    if (this->molecule != 0) {
        THROW_ERROR("Molecule pointer not released");
    }
}

void MolecularGridGenerator::generateGrid(MRGrid<3> &outGrid, Molecule &mol) { 
    this->molecule = &mol;
    GridGenerator<3>::generateGrid(outGrid);
    this->molecule = 0;
}

bool MolecularGridGenerator::splitCheck(GridNode<3> *node) {
    NOT_IMPLEMENTED_ABORT
}


template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
