#include "NuclearPotential.h"
#include "NuclearFunction.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

NuclearPotential::NuclearPotential(double prec, Nuclei &nucs)
        : project(-1.0, MRA->getMaxScale()),
          nuc_func(nucs, prec) {
}

void NuclearPotential::setup(double prec) {
    Timer timer;
    Potential::setup(prec);
    this->project.setPrecision(prec);
    if (this->real == 0) {
        this->real = new FunctionTree<3>(*MRA);
        this->project(*this->real, this->nuc_func);
        this->imag = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
    timer.stop();
    int n = getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Nuclear potential", n, t);
}

void NuclearPotential::clear() {
    this->project.setPrecision(-1.0);
    Potential::clear();
}
