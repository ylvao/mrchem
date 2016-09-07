#include "Density.h"
#include "FunctionTree.h"
#include "TelePrompter.h"

using namespace std;

Density::Density(bool s)
    : spin(s),
      total(0),
      alpha(0),
      beta(0) {
}

Density::Density(const Density &rho)
    : spin(rho.spin),
      total(0),
      alpha(0),
      beta(0) {
}

Density::~Density() {
    if (this->total != 0) MSG_ERROR("Density not properly deallocated");
    if (this->alpha != 0) MSG_ERROR("Density not properly deallocated");
    if (this->beta != 0) MSG_ERROR("Density not properly deallocated");
}

void Density::clear() {
    if (this->total != 0) delete this->total;
    if (this->alpha != 0) delete this->alpha;
    if (this->beta != 0) delete this->beta;
    this->total = 0;
    this->alpha = 0;
    this->beta = 0;
}

int Density::getNNodes() const {
    int nNodes = 0;
    if (this->total != 0) nNodes += this->total->getNNodes();
    if (this->alpha != 0) nNodes += this->alpha->getNNodes();
    if (this->beta != 0) nNodes += this->beta->getNNodes();
    return nNodes;
}

void Density::setDensity(int s, FunctionTree<3> *rho) {
    if (s == Paired) this->total = rho;
    if (s == Alpha) this->alpha = rho;
    if (s == Beta) this->beta = rho;
}

FunctionTree<3>& Density::getDensity(int s) {
    FunctionTree<3> *rho = 0;
    if (s == Paired) rho = this->total;
    if (s == Alpha) rho = this->alpha;
    if (s == Beta) rho = this->beta;
    if (rho == 0) MSG_FATAL("Uninitialized density");
    return *rho;
}

int Density::printTreeSizes() const {
    int nNodes = 0;
    if (this->total != 0) nNodes += this->total->getNNodes();
    if (this->alpha != 0) nNodes += this->alpha->getNNodes();
    if (this->beta != 0) nNodes += this->beta->getNNodes();
    println(0, " Density           " << setw(15) << 1 << setw(25) << nNodes);
    return nNodes;
}

//void Density::setup(double prec) {
//    NOT_IMPLEMENTED_ABORT;
//    Timer timer;

//    FunctionTreeVector<3> sq_vec;
//    FunctionTreeVector<3> sum_vec;
//    this->mult.setPrecision(prec);
//    for (int i = 0; i < this->orbitals->size(); i++) {
//        Orbital &phi_i = this->orbitals->getOrbital(i);
//        double occ = (double) phi_i.getOccupancy();
//        if (phi_i.hasReal()) {
//            sq_vec.push_back(phi_i.re());
//            sq_vec.push_back(phi_i.re());
//            FunctionTree<3> *real_2 = this->mult(sq_vec);
//            sq_vec.clear();
//            sum_vec.push_back(occ, *real_2);
//        }
//        if (phi_i.hasImag()) {
//            sq_vec.push_back(phi_i.im());
//            sq_vec.push_back(phi_i.im());
//            FunctionTree<3> *imag_2 = this->mult(sq_vec);
//            sq_vec.clear();
//            sum_vec.push_back(occ, *imag_2);
//        }
//    }
//    if (this->total == 0) {
//        this->add.setPrecision(prec);
//        this->total = this->add(sum_vec);
//    } else {
//        this->add.setPrecision(-1.0);
//        this->clean.setPrecision(prec);
//        int nNodes = this->clean(*this->total);
//        this->add(*this->total, sum_vec);
//    }

//    for (int i = 0; i < sum_vec.size(); i++) {
//        delete sum_vec[i];
//    }
//    sum_vec.clear();

//    timer.stop();
//    int n = this->total->getNNodes();
//    double t = timer.getWallTime();
//    TelePrompter::printTree(1, "Electron density", n, t);
//    TelePrompter::printDouble(1, "Charge integral", this->total->integrate());
//}
