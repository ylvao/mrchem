#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMMomentum final : public QMOperator {
public:
    QMMomentum(int d, std::shared_ptr<mrcpp::DerivativeOperator<3>> D);

protected:
    const int apply_dir;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
};

class MomentumOperator final : public RankOneTensorOperator<3> {
public:
    MomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
            : p_x(new QMMomentum(0, D))
            , p_y(new QMMomentum(1, D))
            , p_z(new QMMomentum(2, D)) {
        RankOneTensorOperator<3> &p = (*this);
        p[0] = p_x;
        p[1] = p_y;
        p[2] = p_z;
    }

protected:
    std::shared_ptr<QMMomentum> p_x;
    std::shared_ptr<QMMomentum> p_y;
    std::shared_ptr<QMMomentum> p_z;
};

} // namespace mrchem
