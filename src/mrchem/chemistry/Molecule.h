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
class QuadrupoleMoment;
class Polarizability;
class OpticalRotation;
class Magnetizability;
class NMRShielding;
class SpinSpinCoupling;

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
    void initQuadrupoleMoment(const double *o = 0);
    void initPolarizability(double omega = 0, const double *o = 0, bool v = false);
    void initOpticalRotation(double omega = 0, const double *o = 0, bool v = false);
    void initMagnetizability(double omega = 0, const double *o = 0);
    void initNMRShielding(int k, const double *o = 0);
    void initSpinSpinCoupling(int k, int l);

    DipoleMoment &getDipoleMoment();
    QuadrupoleMoment &getQuadrupoleMoment();
    Polarizability &getPolarizability(double omega = 0);
    OpticalRotation &getOpticalRotation(double omega = 0);
    Magnetizability &getMagnetizability(double omega = 0);
    NMRShielding &getNMRShielding(int k);
    SpinSpinCoupling &getSpinSpinCoupling(int k, int l);

protected:
    int charge;
    Nuclei nuclei;

    // Properties
    double COM[3];
    SCFEnergy energy;
    DipoleMoment *dipole;
    QuadrupoleMoment *quadrupole;
    NMRShielding **nmrShielding;
    SpinSpinCoupling ***spinSpinCoupling;
    std::vector<Polarizability *> polarizability;
    std::vector<OpticalRotation *> opticalRotation;
    std::vector<Magnetizability *> magnetizability;

    void allocProperties();
    void freeProperties();

    void clearDipoleMoment();
    void clearQuadrupoleMoment();
    void clearPolarizability();
    void clearOpticalRotation();
    void clearMagnetizability();
    void clearNMRShielding(int k);
    void clearSpinSpinCoupling(int k, int l);

    void calcCenterOfMass();

    void readCoordinateFile(const std::string &file);
    void readCoordinateString(const std::vector<std::string> &coord_str);
};

#endif // MOLECULE_H
