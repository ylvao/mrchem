#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "MomentumOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMMomentum::QMMomentum(int d, std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
        : QMOperator()
        , apply_dir(d)
        , derivative(D) {}

Orbital QMMomentum::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");

    auto dir = this->apply_dir;
    auto &D = *this->derivative;

    Orbital out = inp.paramCopy();

    Timer timer;
    // Calc real part
    if (inp.hasImag()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(out.real(), D, inp.imag(), dir);
        if (inp.conjugate()) out.real().rescale(-1.0);
    }
    // Calc imag part
    if (inp.hasReal()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(out.imag(), D, inp.real(), dir);
        out.imag().rescale(-1.0);
    }
    print_utils::qmfunction(1, "Applied QMMomentum", out, timer);

    return out;
}

Orbital QMMomentum::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
