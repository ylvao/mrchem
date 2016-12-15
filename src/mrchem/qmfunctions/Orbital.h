#ifndef ORBITAL_H
#define ORBITAL_H

#include <complex>

#include "constants.h"

#include "ComplexFunction.h"

class Orbital : public ComplexFunction<3> {
public:
    Orbital(int occ, int s);
    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb);
    virtual ~Orbital() { clear(); }
    void clear(bool free = true);

    int getSpin() const { return this->spin; }
    int getOccupancy() const { return this->occupancy; }
    double getError() const { return this->error; }

    void setSpin(int s) { this->spin = s; }
    void setOccupancy(int occ) { this->occupancy = occ; }
    void setError(double err) { this->error = err; }

    bool isConverged(double prec) const;

    void compare(const Orbital &orb) const;
    int compareSpin(const Orbital &orb) const;
    int compareOccupancy(const Orbital &orb) const;

    std::complex<double> dot(Orbital &ket);
    double getExchangeFactor(const Orbital &orb) const;

    char printSpin() const {
        char sp = 'u';
        if (this->spin == Alpha) sp = 'a';
        if (this->spin == Beta) sp = 'b';
        return sp;
    }

    void send_Orbital(int dest, int tag);
    void Rcv_Orbital(int source, int tag);
    void Isend_Orbital(int dest, int tag);
    void IRcv_Orbital(int source, int tag);

    friend std::ostream& operator<<(std::ostream &o, Orbital &orb) {
        o << std::setw(25) << orb.getSquareNorm();
        o << std::setw(3) << orb.getOccupancy();
        o << std::setw(4) << orb.printSpin();
        o << std::setw(24) << orb.getError() << std::endl;
        return o;
    }

protected:
    int spin;
    int occupancy;
    double error;
};

#endif // ORBITAL_H
