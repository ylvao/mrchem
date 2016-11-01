#ifndef DENSITY_H
#define DENSITY_H

#include "constants.h"

#include "TelePrompter.h"

template<int D> class FunctionTree;

class Density {
public:
    Density(bool s);
    Density(const Density &rho);
    Density &operator=(const Density &rho) { NOT_IMPLEMENTED_ABORT; }
    virtual ~Density();

    bool isSpinDensity() const { return this->spin; }
    int getNNodes() const;
    void clear();

    void setDensity(int s, FunctionTree<3> *rho);
    FunctionTree<3> &getDensity(int s = Paired);

    int printTreeSizes() const;
    void send_Density(int dest, int tag);
    void Rcv_Density(int source, int tag);

    friend class DensityProjector;

protected:
    const bool spin;
    FunctionTree<3> *total;
    FunctionTree<3> *alpha;
    FunctionTree<3> *beta;
};

#endif // DENSITY_H

