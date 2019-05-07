#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NablaOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMNabla::QMNabla(int d, std::shared_ptr<DerivativeOperator<3>> D)
        : QMOperator()
        , apply_dir(d)
        , derivative(D) {}

Orbital QMNabla::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");

    auto dir = this->apply_dir;
    auto &D = *this->derivative;

    Orbital out = inp.paramCopy();

    Timer timer;
    // Calc real part
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(out.real(), D, inp.real(), dir);
    }
    // Calc imag part
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(out.imag(), D, inp.imag(), dir);
        if (inp.conjugate()) out.imag().rescale(-1.0);
    }
    print_utils::qmfunction(1, "Applied QMNabla", out, timer);

    return out;
}

Orbital QMNabla::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
