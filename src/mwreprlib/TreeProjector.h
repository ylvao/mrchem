#ifndef TREEPROJECTOR_H
#define TREEPROJECTOR_H

#include "mwrepr_declarations.h"
#include "TelePrompter.h"

template<int D>
class TreeProjector {
public:
    TreeProjector() { }
    virtual ~TreeProjector() { }

    void calcNodeVector(MRNodeVector &nodeVec) const {
        int nNodes = nodeVec.size();
        printout(10, std::setw(6) << nNodes << " nodes" << std::endl);
#pragma omp parallel shared(nodeVec) firstprivate(nNodes)
{
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = static_cast<MWNode<D> &>(*nodeVec[n]);
            calcNode(node);
        }
}
    }

protected:
    virtual void calcNode(MWNode<D> &node) const = 0;
};

#endif // TREEPROJECTOR_H
