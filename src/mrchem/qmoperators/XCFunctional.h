#ifndef XCFUNCTIONAL_H
#define XCFUNCTIONAL_H

#include <vector>
#include <Eigen/Core>

#include "xcfun.h"

class XCFunctional {
public:
    XCFunctional(bool s, int k = 2);
    virtual ~XCFunctional();

    void evaluate(Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;

    void setFunctional(const std::string &name, double coef = 1.0);
    void setDensityCutoff(double cut) { this->cutoff = cut; }

    int getInputLength() const { return this->inputLength; }
    int getOutputLength() const { return this->outputLength; }

    bool isLDA() const;
    bool isGGA() const;
    bool isSpinSeparated() const { return this->spin; }

private:
    bool spin;
    int order;
    int type;
    int inputLength;
    int outputLength;
    double cutoff;

    xc_functional func;

    void setup();
};

#endif // XCFUNCTIONAL_H
