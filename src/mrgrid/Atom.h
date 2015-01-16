/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#ifndef ATOM_H_
#define ATOM_H_

#include "AtomicElement.h"
#include "TelePrompter.h"

class Atom {
public:
    Atom(const AtomicElement &elm, const double *r = 0) {
        this->element = &elm;
        setCoord(r);
    }

    Atom(const Atom &atom) {
        this->element = &atom.getAtomicElement();
        setCoord(atom.getCoord());
    }

    void setCoord(const double *r) {
        for (int d = 0; d < 3; d++) {
            if (r == 0) {
                this->coord[d] = 0.0;
            } else {
                this->coord[d] = r[d];
            }
        }
    }

    const AtomicElement &getAtomicElement() const { return *element; }
    double getNuclearCharge() const { return this->element->getZ(); }
    const double *getCoord() const { return coord; }

    friend std::ostream& operator<<(std::ostream &o, const Atom &a) {
        o << a.getAtomicElement().getSymbol() << "  [ ";
        for (int i = 0; i < 3; i++) {
            o << std::setw(25) << a.getCoord()[i] << " ";
        }
        o << " ]" << std::endl;
        return o;
    }
private:
    const AtomicElement *element;
    double coord[3];
};

#endif // ATOM_H_
