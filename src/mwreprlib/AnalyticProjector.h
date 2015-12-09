#ifndef ANALYTICPROJECTOR_H
#define ANALYTICPROJECTOR_H

#include "TreeProjector.h"

template<int D>
class AnalyticProjector : public TreeProjector<D> {
public:
    AnalyticProjector(RepresentableFunction<D> &inp_func) {
        this->func = &inp_func;
    }
    ~AnalyticProjector() {
        this->func = 0;
    }

protected:
    RepresentableFunction<D> *func;

    void calcNode(MWNode<D> &node) const;
};

#endif // ANALYTICPROJECTOR_H
