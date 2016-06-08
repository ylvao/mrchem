#ifndef ORBITAL_H
#define ORBITAL_H

#include "constants.h"

#include "TelePrompter.h"

template<int D> class FunctionTree;

class Orbital {
public:
    Orbital(int occ = 2, int s = Paired)
            : spin(s),
              occupancy(occ),
              real(0),
              imag(0) {
        NOT_IMPLEMENTED_ABORT;
    }
    Orbital(const Orbital &orb)
            : spin(orb.spin),
              occupancy(orb.occupancy),
              real(0),
              imag(0) {
        NOT_IMPLEMENTED_ABORT;
    }
    Orbital &operator =(const Orbital &orb) { NOT_IMPLEMENTED_ABORT; }
    virtual ~Orbital() { NOT_IMPLEMENTED_ABORT; }

    int getSpin() const { return this->spin; }
    int getOccupancy() const { return this->occupancy; }
    double getError() const { return this->error; }

    void setSpin(int s) { this->spin = s; }
    void setOccupancy(int occ) { this->occupancy = occ; }
    void setError(double err) { this->error = err; }

    bool isConverged(double prec) const { NOT_IMPLEMENTED_ABORT; }

    FunctionTree<3> &re() { return *this->real; }
    FunctionTree<3> &im() { return *this->imag; }

    void compare(const Orbital &orb) const { NOT_IMPLEMENTED_ABORT; }
    int compareSpin(const Orbital &orb) const { NOT_IMPLEMENTED_ABORT; }
    int compareOccupancy(const Orbital &orb) const { NOT_IMPLEMENTED_ABORT; }

    double dot(Orbital &orb) { NOT_IMPLEMENTED_ABORT; }
    double getSquareNorm() const { NOT_IMPLEMENTED_ABORT; }

    friend std::ostream& operator<<(std::ostream &o, Orbital &orb) {
        char sp = 'u';
        if (orb.getSpin() == Alpha) sp = 'a';
        if (orb.getSpin() == Beta) sp = 'b';
        o << std::setw(25) << sqrt(orb.getSquareNorm());
        o << std::setw(3) << orb.getOccupancy();
        o << std::setw(4) << sp;
        o << std::setw(24) << orb.getError() << std::endl;
        return o;
    }

protected:
    // Parameters
    int spin;
    int occupancy;

    //Data
    double error;
    FunctionTree<3> *real;
    FunctionTree<3> *imag;

//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int version) {
//        ar.register_type(static_cast<Orbital *>(NULL));
//        ar & spin;
//        ar & occupancy;
//        ar & error;
//        ar & boost::serialization::base_object<FunctionTree<3> >(*this);
//    }
};

#endif // ORBITAL_H
