#ifndef ATOM_H
#define ATOM_H

#include "Element.h"

#include <vector>

#define Nuclei std::vector<Nucleus *>

class Nucleus {
public:
    Nucleus(const Element &elm, const double *r = 0)
            : charge(elm.getZ()),
              element(&elm) {
        setCoord(r);
    }
    Nucleus(const Nucleus &nuc)
            : charge(nuc.charge),
              element(nuc.element) {
        setCoord(nuc.coord);
    }
    const Nucleus &operator=(const Nucleus &nuc) {
        if (this != &nuc) {
            this->charge = nuc.charge;
            this->element = nuc.element;
            setCoord(nuc.coord);
        }
        return *this;
    }

    virtual ~Nucleus() { }

    void setCharge(double z) { this->charge = z; }
    void setCoord(const double *r) {
        for (int d = 0; d < 3; d++) {
            if (r != 0) {
                this->coord[d] = r[d];
            } else {
                this->coord[d] = 0.0;
            }
        }
    }

    double getCharge() const { return this->charge; }
    const double *getCoord() const { return this->coord; }
    const Element &getElement() const { return *this->element; }

    friend std::ostream& operator<<(std::ostream &o, const Nucleus &nuc) {
        o << std::endl << *nuc.element << "   Z = ";
        o << nuc.charge << std::endl << "   ";
        for (int i = 0; i < 3; i++) {
            o << nuc.coord[i] << " ";
        }
        o << std::endl;
        return o;
    }

private:
    double charge;
    double coord[3];
    const Element *element;
};

#endif // NUCLEUS_H

