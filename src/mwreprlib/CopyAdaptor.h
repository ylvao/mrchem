#ifndef COPYADAPTOR_H
#define COPYADAPTOR_H

#include "mwrepr_declarations.h"
#include "TreeAdaptor.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor() { }
    CopyAdaptor(const CopyAdaptor<D> &adap) { }
    virtual ~CopyAdaptor() { }

    virtual TreeAdaptor<D> *copy() const { return new CopyAdaptor<D>(*this); }

protected:
    virtual bool splitNode(MWNode<D> &node) const {
        NOT_IMPLEMENTED_ABORT;
    }
};

#endif // COPYADAPTOR_H
