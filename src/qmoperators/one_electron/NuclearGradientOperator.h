#pragma once

#include "QMPotential.h"
#include "analyticfunctions/NuclearGradientFunction.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class NuclearGradientPotential final : public QMPotential {
public:
    NuclearGradientPotential(int d, const Nucleus &nuc, double c)
            : QMPotential(1)
            , func(d, nuc, c) {}

    void setup(double prec) override;
    void clear() override;

protected:
    NuclearGradientFunction func;
};

class NuclearGradientOperator : public RankOneTensorOperator<3> {
public:
    NuclearGradientOperator(const Nucleus &nuc, double c)
            : x_rm3(new NuclearGradientPotential(0, nuc, c))
            , y_rm3(new NuclearGradientPotential(1, nuc, c))
            , z_rm3(new NuclearGradientPotential(2, nuc, c)) {
        RankOneTensorOperator &v = (*this);
        v[0] = x_rm3;
        v[1] = y_rm3;
        v[2] = z_rm3;
    }

protected:
    std::shared_ptr<NuclearGradientPotential> x_rm3;
    std::shared_ptr<NuclearGradientPotential> y_rm3;
    std::shared_ptr<NuclearGradientPotential> z_rm3;
};

} // namespace mrchem
