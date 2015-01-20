#include <fstream>
#include <string>

#include "Intgrl.h"
#include "Atom.h"
#include "PeriodicTable.h"
#include "TelePrompter.h"
#include "MathUtils.h"

using namespace std;

#define GETLINE(X,S)  if (not getline(X,S)) \
    THROW_ERROR("Unexpected end of file while reading basis sets!");

Intgrl::Intgrl(const string &file) {
    fstream ifs;
    ifs.open(file.c_str());
    if (not ifs) {
        THROW_ERROR("Failed to open basis set file: " << file);
    }
    readIntgrlFile(ifs);
    ifs.close();
}

Intgrl::~Intgrl() {
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (this->atoms[i] != 0) {
            delete atoms[i];
        }
    }
    for (unsigned int i = 0; i < basis.size(); i++) {
        if (this->basis[i] != 0) {
            delete basis[i];
        }
    }
}

void Intgrl::readIntgrlFile(iostream &ifs) {
    int nTypes;
    string line;
    GETLINE(ifs, line); //First line is a comment
    GETLINE(ifs, line); //Second line is number of atom types
    istringstream iss(line);
    iss >> nTypes;
    for (int i = 0; i < nTypes; i++) {
        readAtomBlock(ifs);
    }
}

void Intgrl::readAtomBlock(iostream &ifs) {
    double z;
    int nAtoms;
    int nShell;

    ifs >> z;
    ifs >> nAtoms;
    ifs >> nShell;

    int funcsPerShell[nShell];
    for (int i = 0; i < nShell; i++) {
        ifs >> funcsPerShell[i];
    }

    readAtomData(ifs, nAtoms, z);
    AOBasis bas;
    for (int i = 0; i < nShell; i++) {
        for (int j = 0; j < funcsPerShell[i]; j++) {
            readContractionBlock(ifs, bas, i);
        }
    }
    for (int i = 0; i < nAtoms; i++) {
        AOBasis *aoBas = new AOBasis(bas);
        this->basis.push_back(aoBas);
    }
}

void Intgrl::readAtomData(iostream &ifs, int n_atoms, double z) {
    double coord[3];
    string sym;
    double origin[3] = {0.0, 0.0, 0.0};
    //const double *origin = BoundingBox<3>::getWorldBox().getOrigin();
    for (int j = 0; j < n_atoms; j++) {
        ifs >> sym;
        if (sym.size() > 2) {
            sym = sym.erase(2);
        }
        for (int d = 0; d < 3; d++) {
            ifs >> coord[d];
            coord[d] -= origin[d];
        }
        PeriodicTable pt;
        const AtomicElement &element = pt.getAtomicElement(sym.c_str());

        Atom *atom = new Atom(element, coord);
        atom->setNuclearCharge(z);
        this->atoms.push_back(atom);
    }
}

void Intgrl::readContractionBlock(iostream &ifs, AOBasis &basis, int l) {
    if (l > 2) {
        MSG_FATAL("Only s, p and d orbitas are currently supported");
    }
    int nprim, nctr;
    ifs >> nprim;
    ifs >> nctr;

    int start = basis.size();
    int nbas = 0;
    for (int i = 0; i < nctr; i++) {
        AOContraction ctr(l);
        basis.append(ctr);
        nbas += (l + 1) * (l + 2) / 2;
    }

    double expo;
    double coef;
    for (int k = 0; k < nprim; k++) {
        ifs >> expo;
        for (int i = 0; i < nctr; i++) {
            ifs >> coef;
            if (fabs(coef) > MachineZero) {
                basis.getContraction(start + i).append(expo, coef);
            }
        }
    }
}

GaussExp<3> Intgrl::getAtomBasis(int i, bool norm) const {
    assert(i >= 0 and i < this->atoms.size());
    const double *coord = this->atoms[i]->getCoord();
    if (norm) {
        return this->basis[i]->getNormBasis(coord);
    } else {
        return this->basis[i]->getBasis(coord);
    }
}

GaussExp<3> Intgrl::getMolBasis(bool norm) const {
    GaussExp<3> molexp;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        const double *coord = this->atoms[i]->getCoord();
        if (norm) {
            molexp.append(this->basis[i]->getNormBasis(coord));
        } else {
            molexp.append(this->basis[i]->getBasis(coord));
        }
    }
    return molexp;
}
