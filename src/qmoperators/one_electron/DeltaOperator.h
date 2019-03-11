#pragma once

#include "QMPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta final : public QMPotential {
public:
    QMDelta(const mrcpp::Coord<3> &o, double expo);

protected:
    mrcpp::GaussFunc<3> func;

    void setup(double prec) override;
    void clear() override;
};

class DeltaOperator final : public RankZeroTensorOperator {
public:
    DeltaOperator(const mrcpp::Coord<3> &o, double expo = 1.0e6)
            : delta(new QMDelta(o, expo)) {
        RankZeroTensorOperator &d = (*this);
        d = delta;
    }

protected:
    std::shared_ptr<QMDelta> delta;
};

} // namespace mrchem
