/*
 *
 *
 *  \date July, 2010
 *  \author Luca Frediani <luca.frediani@uit.no> \n
 *          CTCC, University of Troms
 *
 *
 *
 * Interface for the XCFun library
 *
 */

#ifndef XCFUN_H_
#define XCFUN_H_

#include <Eigen/Core>
#include <string>

#include <xcfun.h>

template<int D> class NodeIndex;
template<int D> class FunctionTree;

class XCFunctional {
public:
    XCFunctional(bool spinSep, int k = 2);
    virtual ~XCFunctional();

    void setFunctional(const std::string &funcName, double coef = 1.0);

    void evaluate(FunctionTree<3> **input);
    void clear();

    FunctionTree<3> &getOutputFunction(int i);
    void fetchOutputFunction(int i, FunctionTree<3> **output);
    void fetchOutputFunctions(FunctionTree<3> **output);

    int getInputLength() const { return this->inputLength; }
    int getOutputLength() const { return this->outputLength; }

    int getType() const { return this->type; }
    bool isLDA() const;
    bool isGGA() const;
    bool isSpinSeparated() const { return this->spinSeparated; }

    void printInputSizes();
    void printOutputSizes();

private:
    bool spinSeparated;
    int order;
    int type;

    int outMode;
    int inputLength;
    int outputLength;
    int maxInputLength;
    int maxOutputLength;

    xc_functional func;

    FunctionTree<3> **inputFunctions;
    FunctionTree<3> **outputFunctions;

    Eigen::VectorXd ***inputData;
    Eigen::VectorXd ***outputData;

    Eigen::VectorXd &getInputData(int i);
    Eigen::VectorXd &getOutputData(int i);

    void setup();

    Eigen::VectorXd ***allocLocalData(int nFuncs);
    void deleteLocalData();

    template<class T> T **allocPtrArray(int nFuncs);
    void clearInputFunctions();
    void clearOutputFunctions();

    void calcInputData(const NodeIndex<3> &idx);
    void calcInputFunctions(FunctionTree<3> **dens);

    void calcXCValue(Eigen::VectorXd &densityValues);
    void calcXCValues(Eigen::VectorXd &densityValues);
};

#endif /* XCFUN_H_ */
