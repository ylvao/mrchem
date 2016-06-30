#include <fstream>

#include "Molecule.h"
#include "Nucleus.h"
#include "DipoleMoment.h"

using namespace std;
using namespace Eigen;

Molecule::Molecule(const Nuclei &nucs, int c)
        : charge(c),
          nuclei(nucs),
          dipole(0) {
    calcCenterOfMass();
}

Molecule::Molecule(const string &coord_file, int c)
        : charge(c),
          dipole(0) {
    readCoordinateFile(coord_file);
    calcCenterOfMass();
}

Molecule::Molecule(const vector<string> &coord_str, int c)
        : charge(c),
          dipole(0) {
    readCoordinateString(coord_str);
    calcCenterOfMass();
}

Molecule::~Molecule() {
    this->nuclei.clear();
    clearDipoleMoment();
}

void Molecule::clearDipoleMoment() {
    if (this->dipole != 0) {
        delete this->dipole;
        this->dipole = 0;
    }
}

int Molecule::getNElectrons() const {
    int totZ = 0;
    for (int i = 0; i < getNNuclei(); i++) {
        totZ += getNucleus(i).getElement().getZ();
    }
    return totZ - this->charge;
}

void Molecule::initDipoleMoment(const double *o) {
    if (this->dipole != 0) MSG_WARN("Dipole moment already initialized");
    this->dipole = new DipoleMoment();
}

DipoleMoment& Molecule::getDipoleMoment() {
    if (this->dipole == 0) MSG_ERROR("Uninitialized dipole moment");
    return *this->dipole;
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

    int nNuclei;
    string sym;
    double coord[3];
    ifs >> nNuclei;
    for (int i = 0; i < nNuclei; i++) {
        ifs >> sym;
        ifs >> coord[0];
        ifs >> coord[1];
        ifs >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
    ifs.close();
}

void Molecule::readCoordinateString(const vector<string> &coord_str) {
    int nNuclei = coord_str.size();
    string sym;
    double coord[3];
    for (int i = 0; i < nNuclei; i++) {
        stringstream ss;
        ss.str(coord_str[i]);
        ss >> sym;
        ss >> coord[0];
        ss >> coord[1];
        ss >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
}

void Molecule::printGeometry() const {
    TelePrompter::printHeader(0, "Molecule");
    println(0, " Nr  Element             x             y             z      ");
    TelePrompter::printSeparator(0, '-');
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
        printout(0, setw(21) << coord[0]);
        printout(0, setw(14) << coord[1]);
        printout(0, setw(14) << coord[2] << endl);
    }
    TelePrompter::printSeparator(0, '-');
    printout(0, " Center of mass: ");
    printout(0, setw(14) << this->COM[0]);
    printout(0, setw(14) << this->COM[1]);
    printout(0, setw(14) << this->COM[2] << endl);
    TelePrompter::setPrecision(oldPrec);
    TelePrompter::printSeparator(0, '=', 2);
}

void Molecule::printProperties() const {
    println(0, this->energy);
    if (this->dipole != 0) {
        println(0, *this->dipole);
    }
}
