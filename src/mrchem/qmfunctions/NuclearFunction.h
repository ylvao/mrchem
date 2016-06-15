#ifndef NUCLEARFUNCTION_H
#define NUCLEARFUNCTION_H

#include <vector>

#include "RepresentableFunction.h"

class Nuclei;
class Nucleus;

class NuclearFunction : public RepresentableFunction<3> {
public:
    NuclearFunction();
    NuclearFunction(const Nuclei &nucs, double prec);
    virtual ~NuclearFunction() { }

    void push_back(const Nucleus &nuc, double smooth);

    double evalf(const double *r) const;

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

protected:
    const double const_fac;
    std::vector<double> smoothParam;
    std::vector<double> charges;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
};

#endif // NUCLEARFUNCTION_H
