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
#include <Eigen/Core>
#include "GridGenerator.h"

class Molecule;

class MolecularGridGenerator : public GridGenerator<3> {
public:
    MolecularGridGenerator();
    virtual ~MolecularGridGenerator();

    void setAmplitude(int a) { this->amplitude = a; }
    void setWidth(int w) { this->width = w; }
    void setNuclearDependence(int n) { this->nucDep = n; }

    void generateGrid(MRGrid<3> &grid, Molecule &mol);

protected:
    int amplitude; ///< refinement scale for hydrogen 2^{-depth}
    int width;
    int nucDep;
    std::vector<double *> coords;
    std::vector<double > coefs;
    std::vector<double > exps;
    Molecule *molecule;

    bool splitCheck(const GridNode<3> *node);

    void setupRefinementFunction();
    void clearRefinementFunction();
    double evalRefinementFunction(int i, double *r) const;

    bool atomInsideNodeCheck(const GridNode<3> &node);
    bool atomOutsideNodeCheck(const GridNode<3> &node);

    double calcCoef(int Z) const;
    double calcExp(int Z) const;
private:
    static const double widthFac = 1.0e4;
    static const double depthFac = 1.0e1;
};

#endif /* GRID_GENERATOR_H_ */
