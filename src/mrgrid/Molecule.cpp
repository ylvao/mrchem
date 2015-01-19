#include <fstream>

#include "Molecule.h"
#include "Atom.h"
#include "AtomicElement.h"
#include "PeriodicTable.h"
#include "TelePrompter.h"

using namespace std;

Molecule::Molecule(const std::string &coord_file) {
    readCoordinateFile(coord_file);
}

void Molecule::addAtom(const Atom &atom) {
    Atom *newAtom = new Atom(atom);
    this->atoms.push_back(newAtom);
}

void Molecule::readCoordinateFile(const std::string &coord_file) {
    double origin[3] = {0.0, 0.0, 0.0};

    fstream ifs;
    ifs.open(coord_file.c_str());
    if (not ifs) {
        THROW_ERROR("Failed to open basis set file: " << coord_file);
    }

    PeriodicTable pt;
    int nAtoms;
    string sym;
    double coord[3];
    ifs >> nAtoms;
    for (int i = 0; i < nAtoms; i++) {
    ifs >> sym;
    ifs >> coord[0];
    ifs >> coord[1];
    ifs >> coord[2];
    for (int d = 0; d < 3; d++) {
        coord[d] -= origin[d];
    }
    const AtomicElement &element = pt.getAtomicElement(sym.c_str());
    Atom *atom  = new Atom(element, coord);
    this->atoms.push_back(atom);
    }
    ifs.close();
}

Molecule::~Molecule() {
    for (int i = 0; i < getNAtoms(); i++) {
    if (this->atoms[i] != 0) {
        delete this->atoms[i];
        this->atoms[i] = 0;
    }
    }
    this->atoms.clear();
}

void Molecule::print() {
    double origin[3] = {0.0, 0.0, 0.0};

    printout(0, "\n\n========================= Molecule ");
    printout(0, "=========================\n\n");

    int printPrec;
    GET_PRINT_PRECISION(printPrec);
    SET_PRINT_PRECISION(5);

    println(0, " Nr  Element              x             y             z     ");
    println(0, "------------------------------------------------------------");

    int nAtoms = getNAtoms();
    for (int i = 0; i < nAtoms; i++) {
        const Atom &atom = getAtom(i);
        const double *coord = atom.getCoord();
        stringstream symbol;
        symbol << atom.getAtomicElement().getSymbol();
        symbol << "  ";
        printout(0, setw(3) << i+1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, setw(22) << coord[0] + origin[0]);
        printout(0, setw(14) << coord[1] + origin[1]);
        printout(0, setw(14) << coord[2] + origin[2] << endl);
    }
    SET_PRINT_PRECISION(printPrec);
    printout(0, "\n=================================");
    printout(0, "===========================\n\n");
}
