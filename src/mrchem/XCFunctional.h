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

    int getInputLength() const { return this->inputLength; }
    int getOutputLength() const { return this->outputLength; }

    int getType() const { return this->type; }
    bool isLDA() const;
    bool isGGA() const;
    bool isSpinSeparated() const { return this->spinSeparated; }

    void setInputData(int i, Eigen::VectorXd &inpData);
    void calcOutputData(int i, Eigen::VectorXd &outData);

private:
    bool spinSeparated;
    int order;
    int type;

    int inputLength;
    int outputLength;
    int maxInputLength;
    int maxOutputLength;

    xc_functional func;
    Eigen::VectorXd ***inputData;

    void setup();

    Eigen::VectorXd ***allocInputData(int nFuncs);
    void deleteInputData();

    Eigen::VectorXd &getInputData(int i);
};

#endif /* XCFUN_H_ */
