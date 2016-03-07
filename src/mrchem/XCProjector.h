#ifndef XCPROJECTOR_H
#define XCPROJECTOR_H

#include "MWProjector.h"
#include "mwrepr_declarations.h"

class XCFunctional;

class XCProjector : public MWProjector<3> {
public:
    XCProjector(XCFunctional &xc_fun, int k, const MWAdaptor<3> &a = MWAdaptor<3>());
    virtual ~XCProjector();

    void operator()(FunctionTree<3> &out, FunctionTree<3> &inp);

protected:
    int order;
    FunctionTree<3> *inpTree;
    XCFunctional *xcFun;

    void calcNode(MWNode<3> &node);
};

#endif // XCPROJECTOR_H
