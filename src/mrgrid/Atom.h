#ifndef ATOM_H
#define ATOM_H

#include <iostream>

#include "AtomicElement.h"

class Atom {
public:
    Atom(const AtomicElement &elm, const double *r = 0);
    Atom(const Atom &atom);
    const Atom &operator=(const Atom &atom);

    void setCoord(const double *coord);
    void setNuclearCharge(double z) { this->nucCharge = z; }

    const double *getCoord() const { return this->coord; }
    double getNuclearCharge() const { return this->nucCharge; }
    const AtomicElement &getAtomicElement() const { return *this->element; }

    friend std::ostream& operator<<(std::ostream &o, const Atom &a) {
	o << std::endl << a.element->getSymbol()  << "   ";
	for (int i = 0; i < 3; i++) {
	    o << a.coord[i] << " ";
	}
	o << std::endl;
	return o;
    }
private:
    const AtomicElement *element;
    double nucCharge;
    double coord[3];
};

#endif // ATOM_H

