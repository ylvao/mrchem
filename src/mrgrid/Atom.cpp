#include "Atom.h"

Atom::Atom(const AtomicElement &elm, const double *r) {
    this->element = &elm;
    this->nucCharge = this->element->getZ();
    setCoord(r);	
}

Atom::Atom(const Atom &atom) {
    this->element = atom.element;
    setCoord(atom.coord);
    //basis = atom.basis;
}

const Atom &Atom::operator=(const Atom &atom) {
    if (this == &atom) {
	return *this;
    }
    this->element = atom.element;
    this->nucCharge = atom.nucCharge;
    setCoord(atom.coord);
    //basis = atom.basis;
    return *this;
}

void Atom::setCoord(const double *r) {
    for (int d = 0; d < 3; d++) {
	if (r != 0) {
	    this->coord[d] = r[d];
	} else {
	    this->coord[d] = 0.0;
	}
    }
}

//void Atom::setBasis(const AOBasis *abas) {
    //basis = abas;
//}
