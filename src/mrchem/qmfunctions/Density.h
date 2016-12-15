#ifndef DENSITY_H
#define DENSITY_H

#include "constants.h"

template<int D> class FunctionTree;

class Density {
public:
    Density(bool s = false);
    Density(const Density &rho);
    virtual ~Density();
    void clear();

    void setIsSpinDensity(bool s) { this->spin = s; }
    bool isSpinDensity() const { return this->spin; }
    int getNNodes(int type = Paired) const;

    void setDensity(int s, FunctionTree<3> *rho);

    void allocTotal();
    void allocAlpha();
    void allocBeta();

    bool hasTotal() const { if (this->dens_t == 0) return false; return true; }
    bool hasAlpha() const { if (this->dens_a == 0) return false; return true; }
    bool hasBeta() const { if (this->dens_b == 0) return false; return true; }

    FunctionTree<3> &total() { return *this->dens_t; }
    FunctionTree<3> &alpha() { return *this->dens_a; }
    FunctionTree<3> &beta() { return *this->dens_b; }

    const FunctionTree<3> &total() const { return *this->dens_t; }
    const FunctionTree<3> &alpha() const { return *this->dens_a; }
    const FunctionTree<3> &beta() const { return *this->dens_b; }

    void send_Density(int dest, int tag);
    void Rcv_Density(int source, int tag);

protected:
    bool spin;
    FunctionTree<3> *dens_t;
    FunctionTree<3> *dens_a;
    FunctionTree<3> *dens_b;
};

#endif // DENSITY_H

