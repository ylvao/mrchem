#include "DensExp.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "AOBasis.h"
#include "Atom.h"
#include "Intgrl.h"

using namespace std;
using namespace Eigen;

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

DensExp::DensExp(Intgrl &intgrl, Eigen::MatrixXd &D) : densMat(D) {
    this->cartesian = true;
    readAOExpansion(intgrl);
    if (this->size() != this->densMat.rows()) {
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
            prodExp *= this->densMat(i,j);
            densExp.append(prodExp);
        }
    }
    return densExp;
}

double DensExp::evalf(const double *r) const {
    NOT_IMPLEMENTED_ABORT;
    /*
    double val = 0.0;
    for (int i = 0; i < this->size(); i++) {
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

void DensExp::rotate(MatrixXd &U) {
    vector<GaussExp<3> *> tmp;
    int nOrbs = this->orbitals.size();
    for (int i = 0; i < nOrbs; i++) {
        GaussExp<3> *mo = new GaussExp<3>;
    int n = 0;
        for (int j = 0; j < nOrbs; j++) {
            GaussExp<3> tmpExp = *this->orbitals[j];
        //tmpExp.normalize();
            if (fabs(U(i,j)) > MachineZero) {
                tmpExp *= U(i,j);
                mo->append(tmpExp);
        n++;
            //} else {
        //static int nSkip = 0;
        //println(0, "skipping " << nSkip++);
        }
        }
        if (n == 0) {
            MSG_WARN("No contributing orbital");
            GaussExp<3> zeroExp(1);
            GaussFunc<3> zeroFunc(0.0, 0.0);
            zeroExp.setFunc(0, zeroFunc);
            mo->append(zeroExp);
        }
    //mo->normalize();
        tmp.push_back(mo);
    }
    for (int i = 0; i < nOrbs; i++) {
        delete orbitals[i];
        orbitals[i] = tmp[i];
        tmp[i] = 0;
    }
}

void DensExp::readAOExpansion(Intgrl &intgrl) {
    for (int i = 0; i < intgrl.getNAtoms(); i++) {
    Atom &atom = intgrl.getAtom(i);
    AOBasis &aoBasis = intgrl.getAOBasis(i);
    for (int j = 0; j < aoBasis.getNFunc(); j++) {
        GaussExp<3> *ao = new GaussExp<3>(aoBasis.getAO(j, atom.getCoord()));
            this->orbitals.push_back(ao);
        }
    }
    transformToSpherical();
}

void DensExp::transformToSpherical() {
    if (not this->cartesian) {
    return;
    }
    vector<GaussExp<3> *> tmp;
    int nOrbs = this->size();
    int n = 0;
    while (n < nOrbs) {
    int l = getAngularMomentum(n);
    if (l < 2) {
        GaussExp<3> *orb = this->orbitals[n];
        tmp.push_back(orb);
        this->orbitals[n] = 0;
        n++;
    } else if (l == 2) {
        for (int i = 0; i < 5; i++) {
        if (getOrbital(n+i).size() != 1) {
            MSG_FATAL("Cannot handle contracted d orbitals");
        }
        }
        Gaussian<3> &xx = getOrbital(n+0).getFunc(0);
        Gaussian<3> &xy = getOrbital(n+1).getFunc(0);
        Gaussian<3> &xz = getOrbital(n+2).getFunc(0);
        Gaussian<3> &yy = getOrbital(n+3).getFunc(0);
        Gaussian<3> &yz = getOrbital(n+4).getFunc(0);
        Gaussian<3> &zz = getOrbital(n+5).getFunc(0);

        {
        GaussExp<3> *spherical = new GaussExp<3>;
        spherical->append(xy);
        spherical->getFunc(0).setCoef(xy.getCoef());
        spherical->normalize();
        tmp.push_back(spherical);
        }
        {
        GaussExp<3> *spherical = new GaussExp<3>;
        spherical->append(yz);
        spherical->getFunc(0).setCoef(yz.getCoef());
        spherical->normalize();
        tmp.push_back(spherical);
        }
        {
        double coef = 1.0/sqrt(3.0);
        GaussExp<3> *spherical = new GaussExp<3>;
        spherical->append(xx);
        spherical->append(yy);
        spherical->append(zz);
        spherical->getFunc(0).setCoef(-0.5*coef*xx.getCoef());
        spherical->getFunc(1).setCoef(-0.5*coef*yy.getCoef());
        spherical->getFunc(2).setCoef(coef*zz.getCoef());
        spherical->normalize();
        tmp.push_back(spherical);
        }
        {
        GaussExp<3> *spherical = new GaussExp<3>;
        spherical->append(xz);
        spherical->normalize();
        tmp.push_back(spherical);
        }
        {
        GaussExp<3> *spherical = new GaussExp<3>;
        spherical->append(xx);
        spherical->append(yy);
        spherical->getFunc(0).setCoef(0.5*xx.getCoef());
        spherical->getFunc(1).setCoef(-0.5*yy.getCoef());
        spherical->normalize();
        tmp.push_back(spherical);
        }
        n += 6;
    } else {
        MSG_FATAL("Only s, p, and d orbitals are supported");
    }
    }
    for (int i = 0; i < nOrbs; i++) {
        if (this->orbitals[i] != 0) {
        delete this->orbitals[i];
        this->orbitals[i] = 0;
    }
    }
    this->orbitals.clear();
    for (int i = 0; i < tmp.size(); i++) {
        this->orbitals.push_back(tmp[i]);
        tmp[i] = 0;
    }
    this->cartesian = false;
}

int DensExp::getAngularMomentum(int n) const {
    int l = -1;
    GaussExp<3> &gExp = *this->orbitals[n];
    for (int i = 0; i < gExp.size(); i++) {
    const int *pow = gExp.getPower(i);
    int iL = pow[0] + pow[1] + pow[2];
    if (l < 0) {
        l = iL;
    } else if (iL != l) {
        MSG_FATAL("Orbital is not pure angular momentum function");
    }
    }
    return l;
}

/*
bool DensExp::checkSeedNode(MWNode<3> &node) {
    for (unsigned int i = 0; i < this->orbitals.size(); i++) {
        if (getOrbital(i).checkSeedNode(node)) {
            return true;
        }
    }
    return false;
}
*/
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
            for (int j = 0; j < ao.size(); j++) {
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
*/
/*
int DensExp::sortDensTerms(MWNode<3> &node,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &outDens,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &inpVals,
                           Matrix<double, Dynamic, Dynamic, RowMajor> &outVals) {

    Matrix<double, Dynamic, Dynamic, RowMajor> inpDens = this->densMat;
    int n = inpVals.cols();
    int m = 0;
    for (int i = 0; i < this->orbitals.size(); i++) {
        bool contr = false;
        for (int j = 0; j < this->orbitals[i]->size(); j++) {
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
*/
