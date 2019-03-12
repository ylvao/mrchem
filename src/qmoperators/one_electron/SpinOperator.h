#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMSpin final : public QMOperator {
public:
    QMSpin(int d)
            : D(d) {}

protected:
    const int D;

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class SpinOperator final : public RankOneTensorOperator<3> {
public:
    SpinOperator()
            : s_x(0)
            , s_y(1)
            , s_z(2) {
        RankOneTensorOperator &s = *this;
        s[0] = s_x;
        s[1] = s_y;
        s[2] = s_z;
    }

protected:
    QMSpin s_x;
    QMSpin s_y;
    QMSpin s_z;
};

} // namespace mrchem
