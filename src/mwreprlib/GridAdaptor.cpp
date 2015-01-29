#include "GridAdaptor.h"
#include "MRGrid.h"
#include "MRNode.h"
#include "MWNode.h"
#include "NodeIndex.h"
#include "FunctionTree.h"

using namespace std;

template<int D>
GridAdaptor<D>::GridAdaptor(double pr, bool abs) {
    this->prec = pr;
    this->absPrec = abs;
}

template<int D>
void GridAdaptor<D>::adaptGrid(MRGrid<D> &g, FunctionTree<D> &t) {
    this->grid = &g;
    this->tree = &t;

//    println(0, " == Adapting grid");
    MRNodeVector gridVector;
    this->grid->copyEndNodeTable(gridVector);

    NodeIndexSet gridSet;
    for (int i= 0; i < gridVector.size(); i++) {
        const NodeIndex<D> &idx = gridVector[i]->getNodeIndex();
        gridSet.insert(&idx);
    }

    NodeIndexSet splitSet;
    splitNodes(splitSet, gridSet);

    this->grid->clearEndNodeTable();
    this->grid->splitNodes(splitSet);
    this->grid->resetEndNodeTable();

    this->grid = 0;
    this->tree = 0;
}

template<int D>
void GridAdaptor<D>::splitNodes(NodeIndexSet &splitSet, NodeIndexSet &gridSet) {
    typename set<const NodeIndex<D> *>::iterator it;
    for (it = gridSet.begin(); it != gridSet.end(); it++) {
        MRNode<D> *node = this->tree->findNode(**it);
        if (node != 0) {
            MWNode<D> &mwNode = static_cast<MWNode<D> &>(*node);
            if (splitCheck(mwNode)) {
                splitSet.insert(&mwNode.getNodeIndex());
            }
        }
    }
}

template<int D>
bool GridAdaptor<D>::splitCheck(MWNode<D> &node) {
   if (this->prec < 0.0) {
       return false;
   }
   int scale = node.getScale();
   if (scale >= this->tree->getMaxScale()) {
       MSG_INFO("Maximum depth reached: " << scale);
       return false;
   }
   int fact = 1;
   double thr = getWaveletThreshold(fact, scale);
   double w_norm = node.getWaveletNorm();

   if (w_norm > this->prec) {
       return true;
   }
   return false;
}


/** Calculate the threshold for the wavelet norm.
  *
  * Calculates the threshold that has to be met in the wavelet norm in order to
  * guarantee the precision in the function representation. Depends on the
  * square norm of the function and the requested relative accuracy. */
template<int D>
double GridAdaptor<D>::getWaveletThreshold(int factor, int scale) {
    double norm = 1.0;
    if (not this->absPrec) {
        norm = sqrt(this->tree->getSquareNorm());
    }
    double thrs_1 = 2.0 * MachinePrec;
    double thrs_2 = norm * this->prec * pow(2.0, -(0.5 * factor * scale));
    return max(thrs_1, thrs_2);
}


template class GridAdaptor<1>;
template class GridAdaptor<2>;
template class GridAdaptor<3>;
