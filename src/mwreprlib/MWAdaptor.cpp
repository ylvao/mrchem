#include "MWAdaptor.h"
#include "MWNode.h"
#include "constants.h"

using namespace std;

template<int D>
MWAdaptor<D>::MWAdaptor(double pr, bool abs, int scale) {
    this->maxScale = scale;
    this->absPrec = abs;
    this->prec = pr;
}

template<int D>
void MWAdaptor<D>::splitNodeVector(double norm, 
                                    MRNodeVector &nodeVec,
                                    MRNodeVector &splitVec,
                                    MRNodeVector &noSplitVec) {
    for (int n = 0; n < nodeVec.size(); n++) {
        MWNode<D> &node = static_cast<MWNode<D> &>(*nodeVec[n]);
        if (splitCheck(norm, node)) {
            splitVec.push_back(&node);
        } else {
            noSplitVec.push_back(&node);
        }
    }
    nodeVec.clear();
}

template<int D>
bool MWAdaptor<D>::splitCheck(double norm, MWNode<D> &node) {
   if (this->prec < 0.0) {
       return false;
   }
   int scale = node.getScale() + 1;
   if (scale >= this->maxScale) {
       println(10, "Maximum scale reached: " << scale);
       return false;
   }
    if (this->absPrec) {
        norm = 1.0;
    }
   double thrs = getWaveletThreshold(norm, scale);
   double w_norm = node.getWaveletNorm();

   if (w_norm > thrs) {
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
double MWAdaptor<D>::getWaveletThreshold(double norm, int scale) {
    double expo = (0.5 * (scale + 1));
    double thrs_1 = 2.0 * MachinePrec;
    double thrs_2 = norm * this->prec * pow(2.0, -expo);
    return max(thrs_1, thrs_2);
}

template class MWAdaptor<1>;
template class MWAdaptor<2>;
template class MWAdaptor<3>;
