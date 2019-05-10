#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "PositionOperator.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

PositionPotential::PositionPotential(int d, const mrcpp::Coord<3> &o)
        : QMPotential(1) {
    auto f = [d, o](const mrcpp::Coord<3> &r) -> double { return (r[d] - o[d]); };
    this->func.set(f);
}

void PositionPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;
    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
    print_utils::qmfunction(2, "PositionPotential", V, timer);
}

void PositionPotential::clear() {
    free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();    // apply_prec = -1
}

} // namespace mrchem
