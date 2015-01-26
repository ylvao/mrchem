#include "MWProjector.h"
#include "GridAdaptor.h"
#include "MRGrid.h"

template<int D>
MWProjector<D>::MWProjector(MRGrid<D> *startGrid, GridAdaptor<D> *adap) {
    this->grid = startGrid;
    this->adaptor = adap;
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;
