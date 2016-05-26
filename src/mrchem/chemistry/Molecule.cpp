#include <fstream>

#include "Molecule.h"
#include "Nucleus.h"
//#include "DipoleMoment.h"
//#include "QuadrupoleMoment.h"
//#include "Polarizability.h"
//#include "OpticalRotation.h"
//#include "Magnetizability.h"
//#include "NMRShielding.h"
//#include "SpinSpinCoupling.h"
#include "Element.h"
#include "PeriodicTable.h"
//#include "OrbitalSet.h"
//#include "Intgrl.h"
//#include "OrbExp.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

Molecule::Molecule(const Nuclei &nucs, int c)
        : charge(c),
          dipole(0),
          quadrupole(0),
          nmrShielding(0),
          spinSpinCoupling(0) {
    NOT_IMPLEMENTED_ABORT;
    calcCenterOfMass();
    allocProperties();
}

Molecule::Molecule(const string &coord_file, int c)
        : charge(c),
          dipole(0),
          quadrupole(0),
          nmrShielding(0),
          spinSpinCoupling(0) {
    readCoordinateFile(coord_file);
    calcCenterOfMass();
    allocProperties();
}

Molecule::Molecule(const vector<string> &coord_str, int c)
        : charge(c),
          dipole(0),
          quadrupole(0),
          nmrShielding(0),
          spinSpinCoupling(0) {
    readCoordinateString(coord_str);
    calcCenterOfMass();
    allocProperties();
}

Molecule::~Molecule() {
    clearNuclei();
    clearDipoleMoment();
    clearQuadrupoleMoment();
    clearPolarizability();
    clearOpticalRotation();
    clearMagnetizability();
    freeProperties();
}

void Molecule::allocProperties() {
    NOT_IMPLEMENTED_ABORT;
//    int nNucs = this->nuclei.size();
//    this->nmrShielding = new NMRShielding*[nNucs];
//    this->spinSpinCoupling = new SpinSpinCoupling**[nNucs];
//    for (int k = 0; k < nNucs; k++) {
//        this->nmrShielding[k] = 0;
//        this->spinSpinCoupling[k] = new SpinSpinCoupling*[nNucs];
//        for (int l = 0; l < nNucs; l++) {
//            this->spinSpinCoupling[k][l] = 0;
//        }
//    }
}

void Molecule::freeProperties() {
    int nNucs = this->nuclei.size();
    for (int k = 0; k < nNucs; k++) {
        clearNMRShielding(k);
        for (int l = 0; l < nNucs; l++) {
            clearSpinSpinCoupling(k, l);
        }
        delete[] this->spinSpinCoupling[k];
        this->spinSpinCoupling[k] = 0;
    }
    delete[] this->nmrShielding;
    delete[] this->spinSpinCoupling;
    this->nmrShielding = 0;
    this->spinSpinCoupling = 0;
}

void Molecule::clearNuclei() {
    for (int i = 0; i < getNNuclei(); i++) {
        if (this->nuclei[i] != 0) {
            delete this->nuclei[i];
            this->nuclei[i] = 0;
        }
    }
    this->nuclei.clear();
}

void Molecule::clearDipoleMoment() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->dipole != 0) {
//        delete this->dipole;
//        this->dipole = 0;
//    }
}

void Molecule::clearQuadrupoleMoment() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->quadrupole != 0) {
//        delete this->quadrupole;
//        this->quadrupole = 0;
//    }
}

void Molecule::clearNMRShielding(int k) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nmrShielding == 0) MSG_ERROR("Properties not allocated");

//    if (this->nmrShielding[k] != 0) {
//        delete this->nmrShielding[k];
//        this->nmrShielding[k] = 0;
//    }
}

void Molecule::clearSpinSpinCoupling(int k, int l) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->spinSpinCoupling == 0) MSG_ERROR("Properties not allocated");
//    if (this->spinSpinCoupling[k] == 0) MSG_ERROR("Properties not allocated");

//    if (this->spinSpinCoupling[k][l] != 0) {
//        delete this->spinSpinCoupling[k][l];
//        this->spinSpinCoupling[k][l] = 0;
//    }
}

void Molecule::clearPolarizability() {
    NOT_IMPLEMENTED_ABORT;
//    int nPol = this->polarizability.size();
//    for (int i = 0; i < nPol; i++) {
//        if (this->polarizability[i] != 0) {
//            delete this->polarizability[i];
//            this->polarizability[i] = 0;
//        }
//    }
//    this->polarizability.clear();
}

void Molecule::clearOpticalRotation() {
    NOT_IMPLEMENTED_ABORT;
//    int nOpt = this->opticalRotation.size();
//    for (int i = 0; i < nOpt; i++) {
//        if (this->opticalRotation[i] != 0) {
//            delete this->opticalRotation[i];
//            this->opticalRotation[i] = 0;
//        }
//    }
//    this->opticalRotation.clear();
}

void Molecule::clearMagnetizability() {
    NOT_IMPLEMENTED_ABORT;
//    int nMag = this->magnetizability.size();
//    for (int i = 0; i < nMag; i++) {
//        if (this->magnetizability[i] != 0) {
//            delete this->magnetizability[i];
//            this->magnetizability[i] = 0;
//        }
//    }
//    this->magnetizability.clear();
}

int Molecule::getNElectrons() const {
    int totZ = 0;
    for (int i = 0; i < getNNuclei(); i++) {
        totZ += getNucleus(i).getElement().getZ();
    }
    return totZ - this->charge;
}

void Molecule::initDipoleMoment(const double *o) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->dipole != 0) {
//        MSG_WARN("Dipole moment already initialized");
//    }
//    this->dipole = new DipoleMoment(o);
}

void Molecule::initQuadrupoleMoment(const double *o) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->quadrupole != 0) {
//        MSG_WARN("Quadrupole moment already initialized");
//    }
//    this->quadrupole = new QuadrupoleMoment(o);
}

void Molecule::initNMRShielding(int k, const double *o) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nmrShielding == 0) MSG_ERROR("Properties not allocated");

//    const Nucleus &nuc_K = getNucleus(k);
//    if (this->nmrShielding[k] == 0) {
//        this->nmrShielding[k] = new NMRShielding(nuc_K, o);
//    } else {
//        MSG_ERROR("NMR shielding tensor already initialized");
//    }
}

void Molecule::initSpinSpinCoupling(int k, int l) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->spinSpinCoupling == 0) MSG_ERROR("Properties not allocated");
//    if (this->spinSpinCoupling[k] == 0) MSG_ERROR("Properties not allocated");

//    const Nucleus &nuc_K = getNucleus(k);
//    const Nucleus &nuc_L = getNucleus(l);
//    if (this->spinSpinCoupling[k][l] == 0) {
//        this->spinSpinCoupling[k][l] = new SpinSpinCoupling(nuc_K, nuc_L);
//    } else {
//        MSG_ERROR("Spin-spin coupling tensor already initialized");
//    }
}

void Molecule::initPolarizability(double omega, const double *o, bool v) {
    NOT_IMPLEMENTED_ABORT;
//    for (int i = 0; i < this->polarizability.size(); i++) {
//        Polarizability &pol_i = *this->polarizability[i];
//        double w_i = pol_i.getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            MSG_WARN("Polarizability already initialized");
//            return;
//        }
//    }
//    Polarizability *pol = new Polarizability(omega, o, v);
//    this->polarizability.push_back(pol);
}

void Molecule::initOpticalRotation(double omega, const double *o, bool v) {
    NOT_IMPLEMENTED_ABORT;
//    for (int i = 0; i < this->opticalRotation.size(); i++) {
//        OpticalRotation &opt_i = *this->opticalRotation[i];
//        double w_i = opt_i.getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            MSG_WARN("OpticalRotation already initialized");
//            return;
//        }
//    }
//    OpticalRotation *opt = new OpticalRotation(omega, o, v);
//    this->opticalRotation.push_back(opt);
}

void Molecule::initMagnetizability(double omega, const double *o) {
    NOT_IMPLEMENTED_ABORT;
//    for (int i = 0; i < this->magnetizability.size(); i++) {
//        Magnetizability &mag_i = *this->magnetizability[i];
//        double w_i = mag_i.getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            MSG_WARN("Magnetizability already initialized");
//            return;
//        }
//    }
//    Magnetizability *mag = new Magnetizability(omega, o);
//    this->magnetizability.push_back(mag);
}


DipoleMoment& Molecule::getDipoleMoment() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->dipole == 0) {
//        MSG_ERROR("Uninitialized dipole moment");
//    }
//    return *this->dipole;
}

QuadrupoleMoment& Molecule::getQuadrupoleMoment() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->quadrupole == 0) {
//        MSG_ERROR("Uninitialized quadrupole moment");
//    }
//    return *this->quadrupole;
}

Polarizability& Molecule::getPolarizability(double omega) {
    NOT_IMPLEMENTED_ABORT;
//    Polarizability *pol = 0;
//    for (int i = 0; i < this->polarizability.size(); i++) {
//        pol = this->polarizability[i];
//        double w_i = pol->getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            break;
//        } else {
//            pol = 0;
//        }
//    }
//    if (pol == 0) MSG_ERROR("Uninitialized polarizability");
//    return *pol;
}

OpticalRotation& Molecule::getOpticalRotation(double omega) {
    NOT_IMPLEMENTED_ABORT;
//    OpticalRotation *opt = 0;
//    for (int i = 0; i < this->opticalRotation.size(); i++) {
//        opt = this->opticalRotation[i];
//        double w_i = opt->getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            break;
//        } else {
//            opt = 0;
//        }
//    }
//    if (opt == 0) MSG_ERROR("Uninitialized optical rotation");
//    return *opt;
}

Magnetizability& Molecule::getMagnetizability(double omega) {
    NOT_IMPLEMENTED_ABORT;
//    Magnetizability *mag = 0;
//    for (int i = 0; i < this->magnetizability.size(); i++) {
//        mag = this->magnetizability[i];
//        double w_i = mag->getFrequency();
//        if (fabs(w_i - omega) < MachineZero) {
//            break;
//        } else {
//            mag = 0;
//        }
//    }
//    if (mag == 0) MSG_ERROR("Uninitialized magnetizability");
//    return *mag;
}

NMRShielding& Molecule::getNMRShielding(int k) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nmrShielding == 0) MSG_ERROR("Properties not allocated");

//    if (this->nmrShielding[k] == 0) {
//        MSG_ERROR("Uninitialized NMR shielding tensor " << k);
//    }
//    return *this->nmrShielding[k];
}

SpinSpinCoupling& Molecule::getSpinSpinCoupling(int k, int l) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->spinSpinCoupling == 0) MSG_ERROR("Properties not allocated");
//    if (this->spinSpinCoupling[k] == 0) MSG_ERROR("Properties not allocated");

//    if (this->spinSpinCoupling[k][l] == 0) {
//        MSG_ERROR("Uninitialized spin-spin coupling tensor " << k << " " << l);
//    }
//    return *this->spinSpinCoupling[k][l];
}

void Molecule::calcCenterOfMass() {
    this->COM[0] = 0.0;
    this->COM[1] = 0.0;
    this->COM[2] = 0.0;

    double M = 0.0;
    for (int i = 0; i < getNNuclei(); i++) {
        const Nucleus &nuc = getNucleus(i);
        const double *r_i = nuc.getCoord();
        const double m_i = nuc.getElement().getMass();
        for (int d = 0; d < 3; d++) {
            this->COM[d] += r_i[d] * m_i;
        }
        M += m_i;
    }
    for (int d = 0; d < 3; d++) {
        this->COM[d] *= 1.0/M;
    }
}

void Molecule::readCoordinateFile(const string &coord_file) {
    fstream ifs;
    ifs.open(coord_file.c_str());
    if (not ifs) {
        MSG_FATAL("Failed to open basis set file: " << coord_file);
    }

    PeriodicTable pt;
    int nNuclei;
    string sym;
    double coord[3];
    ifs >> nNuclei;
    for (int i = 0; i < nNuclei; i++) {
        ifs >> sym;
        ifs >> coord[0];
        ifs >> coord[1];
        ifs >> coord[2];

        const Element &element = pt.getElement(sym.c_str());
        Nucleus *nuc  = new Nucleus(element, coord);
        this->nuclei.push_back(nuc);
    }
    ifs.close();
}

void Molecule::readCoordinateString(const vector<string> &coord_str) {
    PeriodicTable pt;
    string sym;
    double coord[3];
    for (int i = 0; i < coord_str.size(); i++) {
        stringstream ss;
        ss.str(coord_str[i]);
        ss >> sym;
        ss >> coord[0];
        ss >> coord[1];
        ss >> coord[2];

        const Element &element = pt.getElement(sym.c_str());
        Nucleus *nuc  = new Nucleus(element, coord);
        this->nuclei.push_back(nuc);
    }
}

void Molecule::printGeometry() const {
    println(0, "                                                            ");
    println(0, "========================= Molecule =========================");
    println(0, "                                                            ");
    println(0, " Nr  Element              x             y             z     ");
    println(0, "------------------------------------------------------------");
    println(0, "                                                            ");

    int oldPrec = TelePrompter::setPrecision(5);

    int nNuclei = getNNuclei();
    for (int i = 0; i < nNuclei; i++) {
        const Nucleus &nuc = getNucleus(i);
        const double *coord = nuc.getCoord();
        stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, setw(3) << i+1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, setw(22) << coord[0]);
        printout(0, setw(14) << coord[1]);
        printout(0, setw(14) << coord[2] << endl);
    }
    println(0, "                                                            ");
    println(0, "------------------------------------------------------------");
    println(0, "                                                            ");
    printout(0, " Center Of Mass:  ");
    printout(0, setw(14) << this->COM[0]);
    printout(0, setw(14) << this->COM[1]);
    printout(0, setw(14) << this->COM[2] << endl);
    TelePrompter::setPrecision(oldPrec);
    println(0, "                                                            ");
    println(0, "============================================================");
    println(0, "                                                            ");
    println(0, "                                                            ");
}

void Molecule::printProperties() const {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nmrShielding == 0) MSG_ERROR("Properties not allocated");
//    if (this->spinSpinCoupling == 0) MSG_ERROR("Properties not allocated");

//    println(0, this->energy);
//    if (this->dipole != 0) {
//        println(0, *this->dipole);
//    }
//    if (this->quadrupole != 0) {
//        println(0, *this->quadrupole);
//    }
//    for (int i = 0; i < this->polarizability.size(); i++) {
//        println(0, *this->polarizability[i]);
//    }
//    for (int i = 0; i < this->opticalRotation.size(); i++) {
//        println(0, *this->opticalRotation[i]);
//    }
//    for (int i = 0; i < this->magnetizability.size(); i++) {
//        println(0, *this->magnetizability[i]);
//    }
//    for (int k = 0; k < this->nuclei.size(); k++) {
//        if (this->nmrShielding[k] != 0) {
//            println(0, *this->nmrShielding[k]);
//        }
//        for (int l = 0; l < this->nuclei.size(); l++) {
//            if (this->spinSpinCoupling[k][l] != 0) {
//                println(0, *this->spinSpinCoupling[k][l]);
//            }
//        }
//    }
}
