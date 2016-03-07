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

#include "TelePrompter.h"

class Atom;

class Molecule {
public:
    Molecule() : charge(0), multiplicity(1) { }
    Molecule(const std::string &coord_file);
    ~Molecule();

    void print();
    void addAtom(const Atom &atom);

    void setCharge(int c) { this->charge = c; }
    void setMultiplicity(int m) { this->multiplicity = m; }

    int getCharge() const { return this->charge; }
    int getMultiplicity() const { return this->multiplicity; }

    int getNElectrons() const;
    int getNAtoms() const { return this->atoms.size(); }

    Atom &getAtom(int i) { return *this->atoms[i]; }
    const Atom &getAtom(int i) const { return *this->atoms[i]; }

protected:
    int charge;
    int multiplicity;
    std::vector<Atom *> atoms;

    void readCoordinateFile(const std::string &file);
};

#endif // MOLECULE_H
