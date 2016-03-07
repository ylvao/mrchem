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

#include <vector>
#include "GridGenerator.h"

class Molecule;
class Atom;

class MolecularGridGenerator : public GridGenerator<3> {
public:
    MolecularGridGenerator(double wf = 1.0e2, double df = 1.0e1);
    virtual ~MolecularGridGenerator();

    void setDepth(int d) { if (d < 0) d = 0; this->depth = d; }
    void setWidth(int w) { if (w < 0) w = 0; this->width = w; }
    void setNuclearDependence(int n) { if (n < 0) n = 0; this->nucDep = n; }

    void generateGrid(MRGrid<3> &grid, Molecule &mol);
    void generateGrid(MRGrid<3> &grid, Atom &atom);

protected:
    int depth; ///< refinement scale for hydrogen 2^{-depth}
    int width;
    int nucDep;
    std::vector<double *> coords;
    std::vector<double > coefs;
    std::vector<double > exps;
    Molecule *molecule;

    bool splitCheck(const MRNode<3> *node);

    void setupRefinementFunction();
    void clearRefinementFunction();
    double evalRefinementFunction(int i, double *r) const;

    bool atomInsideNodeCheck(const MRNode<3> &node);
    bool atomOutsideNodeCheck(const MRNode<3> &node);

    double calcCoef(int Z) const;
    double calcExp(int Z) const;
private:
    const double widthFac;
    const double depthFac;
};

#endif /* GRID_GENERATOR_H_ */
