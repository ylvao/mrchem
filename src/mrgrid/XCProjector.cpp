#include "XCProjector.h"
#include "XCFunctional.h"
#include "FunctionTree.h"
#include "MWNode.h"
#include "TelePrompter.h"

using namespace Eigen;

XCProjector::XCProjector(XCFunctional &xc_fun, int k, const MWAdaptor<3> &a) 
     : MWProjector<3>(a, 0) {
    this->order = k;
    this->xcFun = &xc_fun;
}

XCProjector::~XCProjector() {
    this->xcFun = 0;
}

void XCProjector::operator()(FunctionTree<3> &out, FunctionTree<3> &inp) {
    this->inpTree = &inp;
    this->buildTree(out);
    out.mwTransformUp();
    this->inpTree = 0;
}

void XCProjector::calcNode(MWNode<3> &outNode) {
    const NodeIndex<3> &nIdx = outNode.getNodeIndex();
    MWNode<3> &inpNode = static_cast<MWNode<3> &>(this->inpTree->getNode(nIdx));

    inpNode.mwTransform(Reconstruction);
    inpNode.cvTransform(Forward);
    VectorXd inpData = inpNode.getCoefs();
    inpNode.cvTransform(Backward);
    inpNode.mwTransform(Compression);

    VectorXd outData = VectorXd::Zero(inpData.size());
    this->xcFun->setInputData(0, inpData);
    this->xcFun->calcOutputData(this->order, outData);

    outNode.setCoefs(outData);
    outNode.cvTransform(Backward);
    outNode.mwTransform(Compression);
    outNode.setHasCoefs();
    outNode.calcNorms();
}
