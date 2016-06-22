#include "NuclearPotential.h"
#include "NuclearFunction.h"

NuclearPotential::NuclearPotential(double build_prec,
                                   const MultiResolutionAnalysis<3> &mra,
                                   Nuclei &nucs)
        : Potential(mra),
          clean(mra),
          project(mra),
          nuc_func(nucs, build_prec) {
}

NuclearPotential::~NuclearPotential() {
}

void NuclearPotential::setup(double prec) {
    Timer timer;
    timer.restart();
    if (this->real == 0) {
        this->clean.setPrecision(-1.0);
        this->project.setPrecision(prec);
        this->real = this->project(this->nuc_func);
        this->imag = 0;
    } else {
        this->clean.setPrecision(prec);
        this->project.setPrecision(-1.0);
        int nNodes = this->clean(*this->real);
        this->project(*this->real, this->nuc_func);
        this->imag = 0;
    }
    timer.stop();
    int n = getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Nuclear potential", n, t);
    Potential::setup(prec);
}

void NuclearPotential::clear() {
    this->clean.setPrecision(-1.0);
    this->project.setPrecision(-1.0);
    Potential::clear();
}

int NuclearPotential::printTreeSizes() const {
    NOT_IMPLEMENTED_ABORT;
}
