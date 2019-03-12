#pragma once

#include "QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta final : public QMPotential {
public:
    QMDelta(const mrcpp::Coord<3> &o, double expo);

private:
    mrcpp::GaussFunc<3> func;

    void setup(double prec) override;
    void clear() override;
};

class DeltaOperator final : public RankZeroTensorOperator {
public:
    DeltaOperator(const mrcpp::Coord<3> &o, double expo = 1.0e6) {
        delta = std::make_shared<QMDelta>(o, expo);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &d = (*this);
        d = delta;
    }

private:
    std::shared_ptr<QMDelta> delta{nullptr};
};

} // namespace mrchem
