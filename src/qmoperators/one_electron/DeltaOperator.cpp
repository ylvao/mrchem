#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "DeltaOperator.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMDelta::QMDelta(const mrcpp::Coord<3> &o, double expo)
        : QMPotential(1)
        , func(expo, std::pow(expo / MATHCONST::pi, 3.0 / 2.0), o) {}

void QMDelta::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
    print_utils::qmfunction(1, "Delta operator", V, timer);
}

void QMDelta::clear() {
    free(NUMBER::Total);
    clearApplyPrec();
}

} // namespace mrchem
