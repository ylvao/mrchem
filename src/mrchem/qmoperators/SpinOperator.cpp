#include "SpinOperator.h"
#include "MWAdder.h"
#include "GridGenerator.h"
#include "FunctionTreeVector.h"

Orbital* QMSpinX::operator()(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (phi.getSpin() == Paired) NOT_IMPLEMENTED_ABORT;

    MWAdder<3> add(this->apply_prec, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    Timer timer;
    Orbital *Ophi = new Orbital(phi);

    double coef = 0.0;
    if (phi.getSpin() == Alpha) coef = 1.0/2.0;
    if (phi.getSpin() == Beta) coef = 1.0/2.0;

    if (phi.hasReal()) {
        FunctionTreeVector<3> trees;
        trees.push_back(coef, &phi.real());
        Ophi->allocReal();
        grid(Ophi->real(), trees);
        add(Ophi->real(), trees, 0);
    }
    if (phi.hasImag()) {
        FunctionTreeVector<3> trees;
        trees.push_back(coef, &phi.imag());
        Ophi->allocImag();
        grid(Ophi->imag(), trees);
        add(Ophi->imag(), trees, 0);
    }

    if (phi.getSpin() == Alpha) Ophi->setSpin(Beta);
    if (phi.getSpin() == Beta) Ophi->setSpin(Alpha);

    timer.stop();
    int n = Ophi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied spin_x operator", n, t);

    return Ophi;
}

Orbital* QMSpinX::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
}

Orbital* QMSpinY::operator()(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (phi.getSpin() == Paired) NOT_IMPLEMENTED_ABORT;

    MWAdder<3> add(this->apply_prec, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    Timer timer;
    Orbital *Ophi = new Orbital(phi);

    double coef = 0.0;
    if (phi.getSpin() == Alpha) coef = 1.0/2.0;
    if (phi.getSpin() == Beta) coef = -1.0/2.0;

    if (phi.hasReal()) {
        FunctionTreeVector<3> trees;
        trees.push_back(coef, &phi.real());
        Ophi->allocImag();
        grid(Ophi->imag(), trees);
        add(Ophi->imag(), trees, 0);
    }
    if (phi.hasImag()) {
        FunctionTreeVector<3> trees;
        trees.push_back(-coef, &phi.imag());
        Ophi->allocReal();
        grid(Ophi->real(), trees);
        add(Ophi->real(), trees, 0);
    }

    if (phi.getSpin() == Alpha) Ophi->setSpin(Beta);
    if (phi.getSpin() == Beta) Ophi->setSpin(Alpha);

    timer.stop();
    int n = Ophi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied spin_y operator", n, t);

    return Ophi;
}

Orbital* QMSpinY::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
}

Orbital* QMSpinZ::operator()(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (phi.getSpin() == Paired) NOT_IMPLEMENTED_ABORT;

    MWAdder<3> add(this->apply_prec, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    Timer timer;
    Orbital *Ophi = new Orbital(phi);

    double coef = 0.0;
    if (phi.getSpin() == Alpha) coef = 1.0/2.0;
    if (phi.getSpin() == Beta) coef = -1.0/2.0;

    if (phi.hasReal()) {
        FunctionTreeVector<3> trees;
        trees.push_back(coef, &phi.real());
        Ophi->allocReal();
        grid(Ophi->real(), trees);
        add(Ophi->real(), trees, 0);
    }
    if (phi.hasImag()) {
        FunctionTreeVector<3> trees;
        trees.push_back(coef, &phi.imag());
        Ophi->allocImag();
        grid(Ophi->imag(), trees);
        add(Ophi->imag(), trees, 0);
    }

    timer.stop();
    int n = Ophi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied spin_z operator", n, t);

    return Ophi;
}

Orbital* QMSpinZ::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
}

