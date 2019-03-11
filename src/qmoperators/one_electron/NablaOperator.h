#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMNabla final : public QMOperator {
public:
    QMNabla(int d, std::shared_ptr<mrcpp::DerivativeOperator<3>> D);

protected:
    const int apply_dir;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class NablaOperator final : public RankOneTensorOperator<3> {
public:
    NablaOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
            : d_x(new QMNabla(0, D))
            , d_y(new QMNabla(1, D))
            , d_z(new QMNabla(2, D)) {
        RankOneTensorOperator<3> &d = (*this);
        d[0] = d_x;
        d[1] = d_y;
        d[2] = d_z;
    }

protected:
    std::shared_ptr<QMNabla> d_x;
    std::shared_ptr<QMNabla> d_y;
    std::shared_ptr<QMNabla> d_z;
};

} // namespace mrchem
