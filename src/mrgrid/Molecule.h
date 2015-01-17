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

class Atom;

class Molecule {
public:
    Molecule() { }
    Molecule(const std::string &coord_file);
    ~Molecule();

    void print();
    void addAtom(const Atom &atom);

    int getNAtoms() const { return this->atoms.size(); }
    Atom &getAtom(int i) { return *this->atoms[i]; }
    const Atom &getAtom(int i) const { return *this->atoms[i]; }

protected:
    std::vector<Atom *> atoms;

    void readCoordinateFile(const std::string &file);
};

#endif // MOLECULE_H
