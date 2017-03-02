#ifndef COMPLEX_FUNCTION_H
#define COMPLEX_FUNCTION_H

#include <complex>

#include "constants.h"
#include "FunctionTree.h"

template<int D>
class ComplexFunction {
public:
    ComplexFunction(FunctionTree<D> *r = 0, FunctionTree<D> *i = 0) : re(r), im(i) { }
    ComplexFunction(const ComplexFunction &func) { NOT_IMPLEMENTED_ABORT; }
    ComplexFunction<D> &operator=(const ComplexFunction<D> &func);
    virtual ~ComplexFunction() { clearReal(true); clearImag(true); }

    bool hasReal() const { if (this->re == 0) return false; return true; }
    bool hasImag() const { if (this->im == 0) return false; return true; }

    void allocReal();
    void allocImag();

    void clearReal(bool free);
    void clearImag(bool free);

    void setReal(FunctionTree<D> *real) { this->re = real; }
    void setImag(FunctionTree<D> *imag) { this->im = imag; }

    FunctionTree<D> &real() { return *this->re; }
    FunctionTree<D> &imag() { return *this->im; }

    const FunctionTree<D> &real() const { return *this->re; }
    const FunctionTree<D> &imag() const { return *this->im; }

    int getNNodes(int type = Total) const;
    double getSquareNorm(int type = Total) const;
    std::complex<double> dot(ComplexFunction<D> &ket);

    void normalize();
    void operator*=(double c);

    enum NumType { Total, Real, Imag };

protected:
    FunctionTree<D> *re;
    FunctionTree<D> *im;
};

#endif // COMPLEX_FUNCTION_H
