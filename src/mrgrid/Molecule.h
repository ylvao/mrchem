/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>

#include "Atom.h"

class Molecule {
public:
    Molecule() { }
    virtual ~Molecule() {
        int nAtoms = getNAtoms();
	for (int i = 0; i < nAtoms; i++) {
	    if (this->atoms[i] != 0) {
		delete this->atoms[i];
	    }
	    this->atoms[i] = 0;
	}
    }

    int getNAtoms() const { return this->atoms.size(); }
    Atom &getAtom(int i) { return *this->atoms[i]; }
    const Atom &getAtom(int i) const { return *this->atoms[i]; }

    void addAtom(const Atom &atom) {
	Atom *newAtom = new Atom(atom);
	this->atoms.push_back(newAtom);
    }

    friend std::ostream& operator<<(std::ostream &o, const Molecule &m) {
	o << "Molecule contains " << m.getNAtoms() << " atoms: " << std::endl;
	for (int i = 0; i < m.getNAtoms(); i++) {
	    o << std::setw(3) << i << "      " << m.getAtom(i);
	}
	return o;
    }
protected:
    std::vector<Atom *> atoms;
};

#endif // MOLECULE_H_
