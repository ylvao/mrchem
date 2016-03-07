/**
 *
 *
 *  \date July, 2010
 *  \author Luca Frediani <luca.frediani@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 * Interface for the xcfun library of density functionals
 *
 * See xcfun documentation for details
 */

#include "XCFunctional.h"
#include "MWNode.h"
#include "NodeIndex.h"
#include "FunctionTree.h"

using namespace std;
using namespace Eigen;

XCFunctional::XCFunctional(bool spinSep, int k) {
    this->spinSeparated = spinSep;
    this->type = -1.0;
    this->order = k;

    this->inputLength = -1;
    this->outputLength = -1;
    this->maxInputLength = 5;
    this->maxOutputLength = 21;

    this->func = xc_new_functional();
    this->inputData = allocInputData(this->maxInputLength);
}

XCFunctional::~XCFunctional() {
    xc_free_functional(func);
    deleteInputData();
}

VectorXd*** XCFunctional::allocInputData(int nFuncs) {
    if (nFuncs < 0) {
        return 0;
    }
    VectorXd ***data = new VectorXd **[nFuncs];
    int nThreads = omp_get_max_threads();
    for (int i = 0; i < nFuncs; i++) {
        data[i] = new VectorXd *[nThreads];
        for (int j = 0; j < nThreads; j++) {
            data[i][j] = new VectorXd;
        }
    }
    return data;
}

void XCFunctional::deleteInputData() {
    int nThreads = omp_get_max_threads();
    int nInp = this->maxInputLength;
    for (int i = 0; i < nInp; i++) {
        for (int j = 0; j < nThreads; j++) {
            delete this->inputData[i][j];
        }
        delete[] this->inputData[i];
    }
    delete[] this->inputData;
}

void XCFunctional::setInputData(int i, VectorXd &inpData) {
    getInputData(i) = inpData;
}

VectorXd& XCFunctional::getInputData(int i) {
    int thread = omp_get_thread_num();
    return *this->inputData[i][thread];
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

void XCFunctional::setFunctional(const string &funcName, double coef) {
    if (xc_set(this->func, funcName.c_str(), coef)) {
        MSG_ERROR("Invalid functional");
    }

    if (this->isGGA()) {
        return;
    }

    if (funcName == "SLATERx") {
        this->type = XC_LDA;
    } else if (funcName == "VWNc") {
        this->type = XC_LDA;
    } else if (funcName == "VWN5c") {
        this->type = XC_LDA;
    } else if (funcName == "LDA") {
        this->type = XC_LDA;
    } else if (funcName == "TFk") {
        this->type = XC_LDA;
    } else if (funcName == "VWk") {
        this->type = XC_GGA;
    } else if (funcName == "BECKEx") {
        this->type = XC_GGA;
    } else if (funcName == "LYPc") {
        this->type = XC_GGA;
    } else if (funcName == "BLYP") {
        this->type = XC_GGA;
    } else if (funcName == "B3LYP") {
        this->type = XC_GGA;
    } else if (funcName == "PW86c") {
        this->type = XC_GGA;
    } else if (funcName == "PW91x") {
        this->type = XC_GGA;
    } else if (funcName == "PW91k") {
        this->type = XC_GGA;
    } else if (funcName == "PW92c") {
        this->type = XC_GGA;
    } else if (funcName == "BP86") {
        this->type = XC_GGA;
    } else if (funcName == "PBE0") {
        this->type = XC_GGA;
    } else if (funcName == "PBEc") {
        this->type = XC_GGA;
    } else if (funcName == "PBEx") {
        this->type = XC_GGA;
    } else if (funcName == "PBE") {
        this->type = XC_GGA;
    } else {
        MSG_ERROR("Invalid functional");
    }
    setup();
}

void XCFunctional::setup() {
    if (not this->isSpinSeparated()) {
        if (this->isLDA()) {
            xc_eval_setup(func, XC_N, XC_PARTIAL_DERIVATIVES, this->order);
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
            xc_eval_setup(func, XC_N_GNN, XC_PARTIAL_DERIVATIVES, this->order);
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
            xc_eval_setup(func, XC_A_B, XC_PARTIAL_DERIVATIVES, this->order);
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
            xc_eval_setup(func, XC_A_B_GAA_GAB_GBB, XC_PARTIAL_DERIVATIVES, this->order);
            this->inputLength = 5;
            if (this->order == 0) {
                this->outputLength = 1;
            } else if (this->order == 1) {
                this->outputLength = 6;
            } else if (this->order == 2) {
                this->outputLength = 21; //Check ordering of output
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

void XCFunctional::calcOutputData(int k, VectorXd &outData) {
    int nPoints = outData.size();
    double *in = new double[this->inputLength];
    double *out = new double[this->outputLength];
    double thrs = 1.0e-10;
    for (int i = 0; i < nPoints; i++) {
        for (int j = 0; j < this->inputLength; j++) {
            in[j] = this->getInputData(j)[i];
            if (in[j] < thrs) {
                in[j] = thrs;
            }
        }
        if (in[0] < thrs) {
            break;
        }
        xc_eval(this->func, in, out);
        outData[i] = out[k];
        if (k == 2) {
            outData[i] *= in[0];
        }
    }
    delete[] in;
    delete[] out;
}

