#ifndef ORBITALADDER_H
#define ORBITALADDER_H

#include <Eigen/Core>

#include "MWAdder.h"

class Orbital;
class OrbitalVector;

class OrbitalAdder : public MWAdder<3> {
public:
    OrbitalAdder(const MultiResolutionAnalysis<3> &mra, double pr = -1.0)
        : MWAdder(mra, pr) { }
    virtual ~OrbitalAdder() { }

    void operator()(Orbital &phi_ab,
                    double a, Orbital &phi_a,
                    double b, Orbital &phi_b);

    void operator()(Orbital &out,
                    std::vector<double> &coefs,
                    std::vector<Orbital *> &orbs);

    void operator()(OrbitalVector &out,
                    double a, OrbitalVector &inp_a,
                    double b, OrbitalVector &inp_b);

    void operator()(Orbital &out,
                    Eigen::VectorXd &c,
                    OrbitalVector &inp);

    void rotate(OrbitalVector &out,
                Eigen::MatrixXd &U,
                OrbitalVector &inp);

    void inPlace(OrbitalVector &out, double c, OrbitalVector &inp);
    void inPlace(Orbital &out, double c, Orbital &inp);

    using MWAdder<3>::operator();
};

#endif // ORBITALADDER_H
