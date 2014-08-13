/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include <Eigen/Core>

#include "MolecularGridGenerator.h"
#include "Molecule.h"
#include "Atom.h"
#include "MRGrid.h"
#include "MRNode.h"

using namespace Eigen;
using namespace std;

MolecularGridGenerator::MolecularGridGenerator(double wf, double df)
        : widthFac(wf), depthFac(df) {
    this->amplitude = -1;
    this->width = -1;
    this->nucDep = -1;
    this->molecule = 0;
}

MolecularGridGenerator::~MolecularGridGenerator() {
    if (this->molecule != 0) {
        THROW_ERROR("Molecule pointer not released");
    }
    if (this->coords.size() != 0) {
        THROW_ERROR("Coordinate vector not cleared");
    }
    if (this->coefs.size() != 0) {
        THROW_ERROR("Coefficient vector not cleared");
    }
    if (this->exps.size() != 0) {
        THROW_ERROR("Exponent vector not cleared");
    }
}

void MolecularGridGenerator::generateGrid(MRGrid<3> &outGrid, Atom &atom) {
    initGrid(outGrid);
    this->molecule = new Molecule;
    this->molecule->addAtom(atom);
    setupRefinementFunction();
    buildGrid();
    clearRefinementFunction();
    delete this->molecule;
    this->molecule = 0;
    clearGrid();
}

void MolecularGridGenerator::generateGrid(MRGrid<3> &outGrid, Molecule &mol) {
    initGrid(outGrid);
    this->molecule = &mol;
    setupRefinementFunction();
    buildGrid();
    clearRefinementFunction();
    this->molecule = 0;
    clearGrid();
}

void MolecularGridGenerator::setupRefinementFunction() {
    int nAtoms = this->molecule->getNAtoms();
    for (int i = 0; i < nAtoms; i++) {
    Atom &atom = this->molecule->getAtom(i);
    const double *coord = atom.getCoord();
    double *r = new double[3];
    for (int d = 0; d < 3; d++) {
        r[d] = coord[d];
    }
    int z = atom.getNuclearCharge();
    double c = calcCoef(z);
    double e = calcExp(z);
        this->coords.push_back(r);
    this->coefs.push_back(c);
    this->exps.push_back(e);
    }
}

void MolecularGridGenerator::clearRefinementFunction() {
    for (int i = 0; i < this->coords.size(); i++) {
    if (this->coords[i] != 0) {
            delete[] this->coords[i];
        this->coords[i] = 0;
    }
    }
    this->coords.clear();
    this->coefs.clear();
    this->exps.clear();
}

double MolecularGridGenerator::evalRefinementFunction(int i, double *r) const {
    double c = this->coefs[i];
    double e = this->exps[i];
    double gauss = 1.0;
    for (int d = 0; d < 3; d++) {
        double r_0 = this->coords[i][d];
        double relPos = (r[d] - r_0);
        gauss *= exp(-e*relPos*relPos);
    }
    int rootScale = getGrid()->getRootScale();
    int uniform = this->uniformScale;
    if (uniform < rootScale) {
    uniform = rootScale;
    }
    return c*gauss + uniform;
}

bool MolecularGridGenerator::splitCheck(const MRNode<3> *node) {
    if (GridGenerator<3>::splitCheck(node)) {
    //return true;
    }
    if (atomInsideNodeCheck(*node)) {
    //println(0, "insideSplit");
    return true;
    }
    if (atomOutsideNodeCheck(*node)) {
    //println(0, "outsideSplit");
    return true;
    }
    return false;
}

bool MolecularGridGenerator::atomInsideNodeCheck(const MRNode<3> &node) {
    int nAtoms = this->molecule->getNAtoms();
    for (int i = 0; i < nAtoms; i++) {
    if (node.hasCoord(this->coords[i])) {
        if (node.getScale() < this->coefs[i]) {
        return true;
        }
    }
    }
    return false;
}

bool MolecularGridGenerator::atomOutsideNodeCheck(const MRNode<3> &node) {
    if (this->width <= 0) {
    return false;
    }
    double c[3];
    node.getCenter(c);

    int nodeScale = node.getScale();
    int nAtoms = this->molecule->getNAtoms();
    for (int i = 0; i < nAtoms; i++) {
        double reqScale = evalRefinementFunction(i, c);
    if (nodeScale < reqScale) {
        return true;
    }
    }
    return false;
}

double MolecularGridGenerator::calcCoef(int Z) const {
    double nucFac = 1.0;
    if (this->nucDep > 0) {
    double lnZ = log2(1.0*Z);
    nucFac = pow(2.0, (this->nucDep*lnZ)/this->depthFac);
    }
    int rootScale = getGrid()->getRootScale();
    int uniform = this->uniformScale;
    if (uniform < rootScale) {
    uniform = rootScale;
    }
    return nucFac*this->amplitude - uniform;
}

double MolecularGridGenerator::calcExp(int Z) const {
    if (this->width <= 0) {
    return 0.0;
    }
    return this->widthFac/(this->width*this->width);
}

template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
