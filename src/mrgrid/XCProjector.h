#ifndef XCPROJECTOR_H
#define XCPROJECTOR_H

#include "MWProjector.h"

class XCFunctional;

class XCProjector : public MWProjector<3> {
public:
    XCProjector(XCFunctional &xc_fun, int k);
    virtual ~XCProjector();

    void operator()(FunctionTree<3> &out, FunctionTree<3> &inp);

protected:
    int order;
    FunctionTree<3> *inpTree;
    XCFunctional *xcFun;

    void calcWaveletCoefs(MWNode<3> &node);
};

#endif // XCPROJECTOR_H
