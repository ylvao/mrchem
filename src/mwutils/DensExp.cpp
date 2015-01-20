#include <boost/timer.hpp>

#include "MathUtils.h"
#include "DensExp.h"
#include "GaussExp.h"
#include "Intgrl.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

DensExp::DensExp(Intgrl &intgrl, Eigen::MatrixXd &D) : OrbExp(intgrl) {
    this->densityMatrix = D;
    if (this->size() != this->densityMatrix.rows()) {
        MSG_FATAL("Size mismatch density matrix");
    }
}

GaussExp<3> DensExp::getAODensExpansion() {
    GaussExp<3> densExp;
    int nBas = this->orbitals.size();
    for (int i = 0; i < nBas; i++) {
        GaussExp<3> &densTerm1 = getOrbital(i);
        for (int j = 0; j < nBas; j++) {
            GaussExp<3> &densTerm2 = getOrbital(j);
            GaussExp<3> prodExp = densTerm1 * densTerm2;
            prodExp *= 2.0*this->densityMatrix(i,j);
            densExp.append(prodExp);
        }
    }
    return densExp;
}

/** Calculate the scaling and wavelet coefs of all the children, and do the
 * outer product to make the nD-scaling coefs. Since a Gaussian expansion
 * is not separable, we have to do the projection term by term. */
/*
void DensExp::calcWaveletCoefs(MWNode<3> &node) {
    static const int tDim = 8;
    MatrixXd &scaling = node.getMWTree().getTmpScalingCoefs();
    VectorXd &tmpvec = node.getMWTree().getTmpScalingVector();
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int inpos = kp1_d - kp1;
    int scale = node.getNodeIndex().scale();

    int nCoefs = tDim * kp1_d;
    int nAO = this->orbitals.size();

    Matrix<double, Dynamic, Dynamic, RowMajor> aoVals = MatrixXd::Zero(nAO, nCoefs);
    vector<MatrixXd *> quadPts;
    node.getChildrenQuadRoots(quadPts);
    for (int i = 0; i < nAO; i++) {
        GaussExp<3> &ao = getOrbital(i);
        for (int child = 0; child < tDim; child++) {
            int l[3];
            node.calcChildTranslation(child, l);
            for (int j = 0; j < ao.getNFuncs(); j++) {
                const Gaussian<3> &gauss = ao.getFunc(j);
                if (gauss.checkScreen(scale + 1, l)) {
                    continue;
                }
                gauss.evalf(*quadPts[child], scaling);
                tmpvec.segment(inpos, kp1) = scaling.col(0);
                MathUtils::tensorExpandCoefs(3, 0, kp1, kp1_d, scaling, tmpvec);
                aoVals.row(i).segment(child * kp1_d, kp1_d) += tmpvec;
            }
        }
    }
    Matrix<double, Dynamic, Dynamic, RowMajor> sortDens;
    Matrix<double, Dynamic, Dynamic, RowMajor> sortVals;

    int nTerms = sortDensTerms(node, sortDens, aoVals, sortVals);
    if (nTerms > 0) {
#ifdef HAVE_BLAS
        Matrix<double, Dynamic, Dynamic, RowMajor> tmpMat(nTerms, nCoefs);
        double *densData = sortDens.data();
        double *valData = sortVals.data();
        double *tmpData = tmpMat.data();

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nTerms, nCoefs, nTerms,
                    1.0, densData, nTerms, valData, nCoefs, 0.0, tmpData, nCoefs);

        node.allocCoefs();
        for (int i = 0; i < nCoefs; i++) {
            node.getCoefs()[i] = cblas_ddot(nTerms, tmpData, nCoefs, valData, nCoefs);
            tmpData++;
            valData++;
        }
#else
        Matrix<double, Dynamic, Dynamic, RowMajor> densVals = sortDens * sortVals;
        densVals.array() *= sortVals.array();

        node.allocCoefs();
        for (int i = 0; i < nCoefs; i++) {
            node.getCoefs()[i] = densVals.col(i).sum();
        }
#endif
    }

    node.cvTransform(MWNode<3>::Backward);
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
    for (int i = 0; i < quadPts.size(); i++) {
        if (quadPts[i] != 0) {
            delete quadPts[i];
        }
    }
}

int DensExp::sortDensTerms(MWNode<3> &node,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &outDens,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &inpVals,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &outVals) {

    Matrix<double, Dynamic, Dynamic, RowMajor> inpDens = *this->densityMatrix;
    int n = inpVals.cols();
    int m = 0;
    for (int i = 0; i < this->orbitals.size(); i++) {
        bool contr = false;
        for (int j = 0; j < this->orbitals[i]->getNFuncs(); j++) {
            Gaussian<3> &gauss = this->orbitals[i]->getFunc(j);
            if (not gauss.checkScreen(node.getScale(), node.getTranslation())) {
                contr = true;
            }
        }
        if (contr == true) {
            MathUtils::swapRows(inpVals, m, i);
            MathUtils::swapRows(inpDens, m, i);
            MathUtils::swapCols(inpDens, m, i);
            m++;
        }
    }

//	printout(0, "Dimensions: [" << aoVals.rows() << ", " << aoVals.cols());
//	println(0, "] -> [" << sortVals.rows() << ", " << sortVals.cols() << "]");

    outDens = inpDens.block(0,0,m,m);
    outVals = inpVals.block(0,0,m,n);

    return m;
}

bool DensExp::checkSeedNode(MWNode<3> &node) {
    for (unsigned int i = 0; i < this->orbitals.size(); i++) {
        if (getOrbital(i).checkSeedNode(node)) {
            return true;
        }
    }
    return false;
}
*/
double DensExp::evalf(const double *r) const {
    NOT_IMPLEMENTED_ABORT;
    /*
    double val = 0.0;
    for (int i = 0; i < this->getNFuncs(); i++) {
        val += this->getFunc(i).evalf(r);
    }
    return val;
    */
}

void DensExp::calcScreening(double nStdDev) {
    for (int i = 0; i < this->orbitals.size(); i++) {
        this->orbitals[i]->calcScreening(nStdDev);
    }
}

void DensExp::setScreen(bool screen) {
    for (int i = 0; i < this->orbitals.size(); i++) {
        this->orbitals[i]->setScreen(screen);
    }
}
