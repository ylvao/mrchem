#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class DistancePotential final : public QMPotential {
public:
    DistancePotential(double pow, const mrcpp::Coord<3> &r_K, double S);

    void setup(double prec) override;
    void clear() override;

private:
    const double power;
    NuclearFunction func;
};

class DistanceOperator final : public RankZeroTensorOperator {
public:
    DistanceOperator(double pow, const mrcpp::Coord<3> &R, double S = 1.0e-7) {
        r_pow = std::make_shared<DistancePotential>(pow, R, S);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_pow;
    }

private:
    std::shared_ptr<DistancePotential> r_pow{nullptr};
};

} // namespace mrchem
