#ifndef ORBITALADDER_H
#define ORBITALADDER_H

#include <Eigen/Core>

#include "MWAdder.h"
#include "GridGenerator.h"

class Orbital;
class OrbitalVector;

class OrbitalAdder {
public:
    OrbitalAdder(double prec = -1.0) : add(prec) { }
    virtual ~OrbitalAdder() { }

    void setPrecision(double prec) { this->add.setPrecision(prec); }

    void operator()(Orbital &phi_ab,
                    double a, Orbital &phi_a,
                    double b, Orbital &phi_b,
                    bool union_grid);

    void operator()(Orbital &out,
                    std::vector<double> &coefs,
                    std::vector<Orbital *> &orbs,
                    bool union_grid);

    void operator()(OrbitalVector &out,
                    double a, OrbitalVector &inp_a,
                    double b, OrbitalVector &inp_b,
                    bool union_grid);

    void operator()(Orbital &out,
                    const Eigen::VectorXd &c,
                    OrbitalVector &inp,
                    bool union_grid);

    void rotate(OrbitalVector &out, const Eigen::MatrixXd &U, OrbitalVector &inp);
    void rotate(OrbitalVector &out, const Eigen::MatrixXd &U);

    void inPlace(OrbitalVector &out, double c, OrbitalVector &inp);
    void inPlace(Orbital &out, double c, Orbital &inp);

protected:
    MWAdder<3> add;
    GridGenerator<3> grid;
};

#endif // ORBITALADDER_H
