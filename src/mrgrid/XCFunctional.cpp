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

    this->outMode = -1;
    this->inputLength = -1;
    this->outputLength = -1;
    this->maxInputLength = 5;
    this->maxOutputLength = 21;

    this->func = xc_new_functional();

    this->inputFunctions = allocPtrArray<FunctionTree<3> >(this->maxInputLength);;
    this->outputFunctions = allocPtrArray<FunctionTree<3> >(this->maxOutputLength);;

    this->inputData = allocLocalData(this->maxInputLength);
    this->outputData = allocLocalData(this->maxOutputLength);
}

XCFunctional::~XCFunctional() {
    xc_free_functional(func);

    clearInputFunctions();
    clearOutputFunctions();

    delete[] this->inputFunctions;
    delete[] this->outputFunctions;

    deleteLocalData();
}

void XCFunctional::clear() {
    this->outMode = -1;
    clearInputFunctions();
    clearOutputFunctions();
}

VectorXd*** XCFunctional::allocLocalData(int nFuncs) {
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

template<class T>
T** XCFunctional::allocPtrArray(int nFuncs) {
    if (nFuncs < 0) {
        return 0;
    }
    T **funcs = new T*[nFuncs];
    for (int i = 0; i < nFuncs; i++) {
        funcs[i] = 0;
    }
    return funcs;
}

void XCFunctional::deleteLocalData() {
    int nThreads = omp_get_max_threads();
    int nInp = this->maxInputLength;
    for (int i = 0; i < nInp; i++) {
        for (int j = 0; j < nThreads; j++) {
            delete this->inputData[i][j];
        }
        delete[] this->inputData[i];
    }
    delete[] this->inputData;

    int nOut = this->maxOutputLength;
    for (int i = 0; i < nOut; i++) {
        for (int j = 0; j < nThreads; j++) {
            delete this->outputData[i][j];
        }
        delete[] this->outputData[i];
    }
    delete[] this->outputData;
}

void XCFunctional::clearInputFunctions() {
    for (int i = 0; i < this->maxInputLength; i++) {
        this->inputFunctions[i] = 0;
    }
}

void XCFunctional::clearOutputFunctions() {
    for (int i = 0; i < this->maxOutputLength; i++) {
        if (this->outputFunctions[i] != 0) {
            delete this->outputFunctions[i];
            this->outputFunctions[i] = 0;
        }
    }
}

VectorXd& XCFunctional::getInputData(int i) {
    int thread = omp_get_thread_num();
    return *this->inputData[i][thread];
}

VectorXd& XCFunctional::getOutputData(int i) {
    int thread = omp_get_thread_num();
    return *this->outputData[i][thread];
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

void XCFunctional::printInputSizes() {
    println(0, endl << "Input functions");
    for (int i = 0; i < this->maxInputLength; i++) {
        if (this->inputFunctions[i] != 0) {
            println(0, i << ": " << this->inputFunctions[i]->getNNodes());
        }
    }
    println(0, "");
}

void XCFunctional::printOutputSizes() {
    println(0, endl << "Output functions");
    for (int i = 0; i < this->maxOutputLength; i++) {
        if (this->outputFunctions[i] != 0) {
            println(0, i << ": " << this->outputFunctions[i]->getNNodes());
        }
    }
    println(0, "");
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

void XCFunctional::evaluate(FunctionTree<3> **input) {
    clearInputFunctions();
    clearOutputFunctions();
    for (int i = 0; i < this->inputLength; i++) {
        if (input[i] == 0) {
            MSG_ERROR("Invalid input");
        }
        this->inputFunctions[i] = input[i];
    }
    for (int i = 0; i < this->outputLength; i++) {
        this->outMode = i;
        this->outputFunctions[i] = new FunctionTree<3>;
        NOT_IMPLEMENTED_ABORT;
//        this->outputFunctions[i]->applyFunctional(*this);
    }
    printInputSizes();
    printOutputSizes();
    clearInputFunctions();
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
FunctionTree<3>& XCFunctional::getOutputFunction(int i) {
    if (i < 0 or i > this->outputLength) {
        MSG_ERROR("Invalid output function");
    }
    if (this->outputFunctions[i] == 0) {
        MSG_ERROR("Output function is not calculated");
    }
    return *this->outputFunctions[i];
}

void XCFunctional::fetchOutputFunction(int i, FunctionTree<3> **output) {
    if (i < 0 or i > this->outputLength) {
        MSG_ERROR("Invalid output function");
    }
    if (this->outputFunctions[i] == 0) {
        MSG_ERROR("Output function is not calculated");
    }
    output[i] = this->outputFunctions[i];
    this->outputFunctions[i] = 0;
}

void XCFunctional::fetchOutputFunctions(FunctionTree<3> **output) {
    for (int i = 0; i < this->outputLength; i++) {
        if (this->outputFunctions[i] == 0) {
            MSG_ERROR("Output function is not calculated");
        }
        output[i] = this->outputFunctions[i];
        this->outputFunctions[i] = 0;
    }
}

//bool XCFunctional::checkSeedNode(MWNode<3> &node) {
//    if (this->inputFunctions[0] == 0) {
//        MSG_ERROR("No density for seeding");
//    }
//    return this->inputFunctions[0]->checkSeedNode(node);
//}

//void XCFunctional::calcWaveletCoefs(MWNode<3> &node) {
//    calcInputData(node.getNodeIndex());
//    calcXCValue(node.getCoefs());
//    node.cvTransform(MWNode<3>::Backward);
//    node.mwTransform(Compression);
//    node.setHasCoefs();
//    node.calcNorms();
//}

void XCFunctional::calcInputData(const NodeIndex<3> &idx) {
    NOT_IMPLEMENTED_ABORT;
//    int nInp = this->inputLength;
//    for (int i = 0; i < nInp; i++) {
//        MRNode<3> &mrNode = this->inputFunctions[i]->getNode(idx);
//        MWNode<3> &mwNode = static_cast<MWNode<3> &>(mrNode);
//        VectorXd &locData = this->getInputData(i);
//        mwNode.mwTransform(Reconstruction);
//        mwNode.cvTransform(Forward);
//        locData = mwNode.getCoefs();
//        mwNode.cvTransform(Backward);
//        mwNode.mwTransform(Compression);
//    }
}

void XCFunctional::calcXCValue(VectorXd &outValues) {
    NOT_IMPLEMENTED_ABORT;
//    int nPoints = outValues.size();
//    double *in = new double[this->inputLength];
//    double *out = new double[this->outputLength];
//    double thrs = 1.0e-12;
//    for (int i = 0; i < nPoints; i++) {
//        for (int j = 0; j < this->inputLength; j++) {
//            in[j] = this->getInputData(j)[i];
//            if (in[j] < 0.0) {
//                in[j] = 0.0;
//            }
//        }
//        if (in[0] < thrs) {
//            outValues = VectorXd::Zero(nPoints);
//            break;
//        }
//        xc_eval(this->func, in, out);
//        outValues[i] = out[this->outMode];
//    }
//    delete[] in;
//    delete[] out;
}

void XCFunctional::calcXCValues(VectorXd &IOValues) {
    NOT_IMPLEMENTED_ABORT
}
