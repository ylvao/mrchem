#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SpinOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

Orbital QMSpin::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    ComplexDouble coef(0.0, 0.0);
    switch (inp.spin()) {
        case SPIN::Alpha:
            if (this->D == 0) coef = ComplexDouble(0.5, 0.0);
            if (this->D == 1) coef = ComplexDouble(0.0, 0.5);
            if (this->D == 2) coef = ComplexDouble(0.5, 0.0);
            break;
        case SPIN::Beta:
            if (this->D == 0) coef = ComplexDouble(0.5, 0.0);
            if (this->D == 1) coef = ComplexDouble(0.0, -0.5);
            if (this->D == 2) coef = ComplexDouble(-0.5, 0.0);
            break;
        default:
            MSG_ABORT("Cannot apply spin operator on paired orbital");
    }

    Orbital out = inp.paramCopy();
    qmfunction::deep_copy(out, inp);
    out.rescale(coef);

    // Flip spin for s_x and s_y
    if (this->D == 0 or this->D == 1) {
        if (inp.spin() == SPIN::Alpha) out.setSpin(SPIN::Beta);
        if (inp.spin() == SPIN::Beta) out.setSpin(SPIN::Alpha);
    }

    return out;
}

Orbital QMSpin::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
