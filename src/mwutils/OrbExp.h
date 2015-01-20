#ifndef ORBEXP_H
#define ORBEXP_H

#include <Eigen/Core>
#include <vector>

class Intgrl;
template<int D> class GaussExp;

class OrbExp {
public:
    OrbExp(Intgrl &intgrl);
    virtual ~OrbExp();

    int size() const { return orbitals.size(); }
    int getAngularMomentum(int n) const;

    GaussExp<3> &getOrbital(int n) { return *orbitals[n]; }

    void rotate(const Eigen::MatrixXd &U);

protected:
    bool cartesian;
    std::vector<GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

#endif // ORBEXP_H
