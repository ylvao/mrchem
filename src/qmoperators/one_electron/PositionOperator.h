#pragma once

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class PositionPotential final : public QMPotential {
public:
    PositionPotential(int d, const mrcpp::Coord<3> &o);

protected:
    mrcpp::AnalyticFunction<3> func;

    void setup(double prec) override;
    void clear() override;
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const mrcpp::Coord<3> &o = {0.0, 0.0, 0.0})
            : r_x(new PositionPotential(0, o))
            , r_y(new PositionPotential(1, o))
            , r_z(new PositionPotential(2, o)) {
        RankOneTensorOperator &r = (*this);
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
    }
    ~PositionOperator() override = default;

protected:
    std::shared_ptr<PositionPotential> r_x;
    std::shared_ptr<PositionPotential> r_y;
    std::shared_ptr<PositionPotential> r_z;
};

} // namespace mrchem
