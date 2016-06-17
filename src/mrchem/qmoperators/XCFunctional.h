#ifndef XCFUNCTIONAL_H
#define XCFUNCTIONAL_H

#include <Eigen/Core>

#include "xcfun.h"

class XCFunctional {
public:
    XCFunctional(bool s);
    virtual ~XCFunctional();

    void setFunctional(const std::string &name, double coef = 1.0);

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength(int k) const { return xc_output_length(this->functional, k); }

    bool isLDA() const { return (xc_get_type(this->functional) == XC_LDA); }
    bool isGGA() const { return (xc_get_type(this->functional) == XC_GGA); }
    bool isSpinSeparated() const { return this->spin; }

    void evaluate(int k, Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;

private:
    bool spin;
    xc_functional functional;

    int getParamFromName(const std::string &name);
};

#endif // XCFUNCTIONAL_H
