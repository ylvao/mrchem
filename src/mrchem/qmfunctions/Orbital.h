#ifndef ORBITAL_H
#define ORBITAL_H

#include <complex>

#include "constants.h"

#include "TelePrompter.h"

template<int D> class FunctionTree;

class Orbital {
public:
    Orbital(int occ, int s);
    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb);
    virtual ~Orbital() { clear(); }
    void clear(bool free = true);

    int getNNodes() const;
    int getSpin() const { return this->spin; }
    int getOccupancy() const { return this->occupancy; }
    double getError() const { return this->error; }

    void setSpin(int s) { this->spin = s; }
    void setOccupancy(int occ) { this->occupancy = occ; }
    void setError(double err) { this->error = err; }

    bool isConverged(double prec) const;

    bool hasReal() const { if (this->real == 0) return false; return true; }
    bool hasImag() const { if (this->imag == 0) return false; return true; }

    FunctionTree<3> &re() { return *this->real; }
    FunctionTree<3> &im() { return *this->imag; }

    void compare(const Orbital &orb) const;
    int compareSpin(const Orbital &orb) const;
    int compareOccupancy(const Orbital &orb) const;

    std::complex<double> dot(Orbital &ket);
    double getSquareNorm() const;
    double getExchangeFactor(const Orbital &orb) const;

    void normalize();
    void orthogonalize(Orbital &orb);

    char printSpin() const {
        char sp = 'u';
        if (this->spin == Alpha) sp = 'a';
        if (this->spin == Beta) sp = 'b';
        return sp;
    }

    friend std::ostream& operator<<(std::ostream &o, Orbital &orb) {
        o << std::setw(25) << orb.getSquareNorm();
        o << std::setw(3) << orb.getOccupancy();
        o << std::setw(4) << orb.printSpin();
        o << std::setw(24) << orb.getError() << std::endl;
        return o;
    }

    friend class MomentumOperator;
    friend class Potential;
    friend class ExchangePotential;
    friend class OrbitalAdder;
    friend class OrbitalMultiplier;
    friend class OrbitalProjector;
    friend class HelmholtzOperatorSet;

protected:
    int spin;
    int occupancy;

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
