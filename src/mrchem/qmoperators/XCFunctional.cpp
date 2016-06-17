#include "XCFunctional.h"
#include "TelePrompter.h"
#include "Timer.h"

#include "constants.h"

using namespace std;
using namespace Eigen;

XCFunctional::XCFunctional(bool s, int k)
        : spin(s),
          order(k),
          type(XC_Undefined),
          inputLength(-1),
          outputLength(-1),
          cutoff(-1.0) {
    this->func = xc_new_functional();
}

XCFunctional::~XCFunctional() {
    xc_free_functional(func);
}

bool XCFunctional::isLDA() const {
    if (this->type == XC_LDA) {
        return true;
    } else {
        return false;
    }
}

bool XCFunctional::isGGA() const {
    if (this->type == XC_GGA) {
        return true;
    } else {
        return false;
    }
}

void XCFunctional::setFunctional(const string &name, double coef) {
    if (xc_set(this->func, name.c_str(), coef)) {
        MSG_ERROR("Invalid functional");
    }

    if (this->isGGA()) {
        return;
    }

    if (name == "SLATERx") {
        this->type = XC_LDA;
    } else if (name == "VWNc") {
        this->type = XC_LDA;
    } else if (name == "VWN5c") {
        this->type = XC_LDA;
    } else if (name == "LDA") {
        this->type = XC_LDA;
    } else if (name == "PW92c") {
        this->type = XC_LDA;
    } else if (name == "PZ81c") {
        this->type = XC_LDA;
    } else if (name == "TFk") {
        this->type = XC_LDA;
    } else if (name == "VWk") {
        this->type = XC_GGA;
    } else if (name == "BECKEx") {
        this->type = XC_GGA;
    } else if (name == "LYPc") {
        this->type = XC_GGA;
    } else if (name == "BLYP") {
        this->type = XC_GGA;
    } else if (name == "B3LYP") {
        this->type = XC_GGA;
    } else if (name == "PW91x") {
        this->type = XC_GGA;
    } else if (name == "PW91c") {
        this->type = XC_GGA;
    } else if (name == "PW91k") {
        this->type = XC_GGA;
    } else if (name == "BP86") {
        this->type = XC_GGA;
    } else if (name == "PBE0") {
        this->type = XC_GGA;
    } else if (name == "PBEc") {
        this->type = XC_GGA;
    } else if (name == "REVPBEc") {
        this->type = XC_GGA;
    } else if (name == "PBEx") {
        this->type = XC_GGA;
    } else if (name == "PBE") {
        this->type = XC_GGA;
    } else {
        MSG_ERROR("Invalid functional");
    }
    setup();
}

void XCFunctional::setup() {
    if (not this->isSpinSeparated()) {
        if (this->isLDA()) {
            xc_eval_setup(this->func,
                          XC_N,
                          XC_PARTIAL_DERIVATIVES,
                          this->order);
            this->inputLength = 1;
            if (this->order == 0) {
                this->outputLength = 1;
            } else if (this->order == 1) {
                this->outputLength = 2;
            } else if (this->order == 2) {
                this->outputLength = 3;
            } else {
                MSG_ERROR("Order > 2 not supported");
            }
        } else if (this->isGGA()) {
            xc_eval_setup(this->func,
                          XC_N_GNN,
                          XC_PARTIAL_DERIVATIVES,
                          this->order);
            this->inputLength = 2;
            if (this->order == 0) {
                this->outputLength = 1;
            } else if (this->order == 1) {
                this->outputLength = 3;
            } else if (this->order == 2) {
                this->outputLength = 6;
            } else {
                MSG_ERROR("Order > 2 not supported");
            }
        } else {
            MSG_ERROR("Invalid functional type");
        }
    } else {
        if (this->isLDA()) {
            xc_eval_setup(this->func,
                          XC_A_B,
                          XC_PARTIAL_DERIVATIVES,
                          this->order);
            this->inputLength = 2;
            if (this->order == 0) {
                this->outputLength = 1;
            } else if (this->order == 1) {
                this->outputLength = 3;
            } else if (this->order == 2) {
                this->outputLength = 6;
            } else {
                MSG_ERROR("Order > 2 not supported");
            }
        } else if (this->isGGA()) {
            xc_eval_setup(this->func,
                          XC_A_B_GAA_GAB_GBB,
                          XC_PARTIAL_DERIVATIVES,
                          this->order);
            this->inputLength = 5;
            if (this->order == 0) {
                this->outputLength = 1;
            } else if (this->order == 1) {
                this->outputLength = 6;
            } else if (this->order == 2) {
                this->outputLength = 21;
            } else {
                MSG_ERROR("Order > 2 not supported");
            }
        } else {
            MSG_ERROR("Invalid functional type");
        }
    }
}

/** Computes the alpha and beta exchange-correlation potentials
 * from the xcfun output functions. For LDA's these are the second
 * and third output functions, respectively. For GGA's the potentials
 * must be computed through
 * \f$ v_{xc}^\sigma = \frac{\partial F_{xc}}{\partial \rho^\sigma(r)}
 *  - \nabla\cdot\frac{\partial F_{xc}}{\partial(\nabla\rho^\sigma)} \f$
 *
 * XCFunctional output:
 *
 * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \rho_x},
 *  \frac{\partial F_{xc}}{\partial \rho_y},
 *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
 *
 * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
 *  \right) \f$
 */

void XCFunctional::evaluate(MatrixXd &inp, MatrixXd &out) const {
    if (inp.cols() != this->inputLength) MSG_ERROR("Invalid input");

    Timer timer;
    timer.restart();
    TelePrompter::printHeader(0, "Evaluating XC functional");

    int nPts = inp.rows();
    out = MatrixXd::Zero(nPts, this->outputLength);

    for (int i = 0; i < nPts; i++) {
        if (inp(i,0) >= this->cutoff) {
            xc_eval(this->func, inp.row(i).data(), out.row(i).data());
        } else {
            out.row(i).setZero();
        }
    }
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}
