/**
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>
#include <Eigen/Core>

#include "SCFEnergy.h"
#include "Nucleus.h"

class DipoleMoment;

class Molecule {
public:
    Molecule(const Nuclei &nucs, int c = 0);
    Molecule(const std::string &coord_file, int c = 0);
    Molecule(const std::vector<std::string> &coord_str, int c = 0);
    virtual ~Molecule();

    int getCharge() const { return this->charge; }
    int getNNuclei() const { return this->nuclei.size(); }
    int getNElectrons() const;

    Nuclei &getNuclei() { return this->nuclei; }
    Nucleus &getNucleus(int i) { return this->nuclei[i]; }
    const Nucleus &getNucleus(int i) const { return this->nuclei[i]; }

    void printGeometry() const;
    void printProperties() const;

    SCFEnergy &getSCFEnergy() { return this->energy; }
    const double *getCenterOfMass() const { return this->COM; }

    void initDipoleMoment(const double *o = 0);

    DipoleMoment &getDipoleMoment();

protected:
    int charge;
    Nuclei nuclei;

    // Properties
    double COM[3];
    SCFEnergy energy;
    DipoleMoment *dipole;

    void calcCenterOfMass();

    void clearDipoleMoment();

    void readCoordinateFile(const std::string &file);
    void readCoordinateString(const std::vector<std::string> &coord_str);
};

#endif // MOLECULE_H
