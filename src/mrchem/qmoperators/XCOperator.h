#ifndef XCOPERATOR_H
#define XCOPERATOR_H

#include <Eigen/Core>

#include "QMOperator.h"
#include "Density.h"
#include "MWAdder.h"
#include "MWMultiplier.h"
#include "DensityProjector.h"
#include "DerivativeOperator.h"
#include "Timer.h"

class XCFunctional;
class OrbitalVector;
class Potential;
template<int D> class FunctionTree;
template<int D> class MultiResolutionAnalysis;

class XCOperator : public QMOperator {
public:
    XCOperator(int k,
               double build_prec,
               const MultiResolutionAnalysis<3> &mra,
               XCFunctional &func,
               OrbitalVector &phi);
    virtual ~XCOperator();

    double getEnergy() const { return this->energy; }

    virtual void setup(double prec);
    virtual void clear();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    const int order;
    XCFunctional *functional;
    MWAdder<3> add;
    MWMultiplier<3> mult;
    DensityProjector project;
    DerivativeOperator<3> derivative;  ///< Derivative operator for GGAs

    Density density_0;              ///< Unperturbed density
    Density **gradient_0;           ///< Unperturbed density gradient
    OrbitalVector *orbitals_0;      ///< Unperturbed orbitals

    double energy;                  ///< XC energy
    Potential *potential[3];        ///< The actual operator [tot, alpha, beta]

    FunctionTree<3> **xcInput;      ///< XCFun input
    FunctionTree<3> **xcOutput;     ///< XCFun output

    void setupXCInput();
    void setupXCOutput();

    void clearXCInput();
    void clearXCOutput();

    void calcDensity();
    Density **calcDensityGradient(Density &rho);

    virtual void calcPotential() = 0;
    bool cropPotential(double prec);

    void calcEnergy();
    void evaluateXCFunctional();

    void compressTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    void compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    FunctionTree<3>* calcDotProduct(FunctionTreeVector<3> &vec_a,
                                    FunctionTreeVector<3> &vec_b);
    FunctionTree<3>* calcGradDotPotDensVec(FunctionTree<3> &pot,
                                           FunctionTreeVector<3> &dens);
//    Potential* calcPotDensVecDotDensVec(Potential *pot, Density **dens_1, Density **dens_2);

    template<class T>
    int sumNodes(T **trees, int nTrees) {
        int nNodes = 0;
        for (int i = 0; i < nTrees; i++) {
            if (trees[i] != 0) {
                nNodes += trees[i]->getNNodes();
            }
        }
        return nNodes;
    }

//    template<class T>
//    T* calcDivergence(T **vec) {
//        Timer timer;
//        if (vec == 0) MSG_ERROR("Input vector not initialized");

//        vector<FunctionTree<3> *> trees;
//        for (int d = 0; d < 3; d++) {
//            timer.resume();
//            if (vec[d] == 0) MSG_ERROR("Invalid input");
//            FunctionTree<3> *tree = new FunctionTree<3>;
//            this->derivative->setApplyDir(d);
//            this->derivative->apply(*tree, *vec[d]);
//            trees.push_back(tree);
//            double time = timer.elapsed();
//            int nNodes = tree->getNNodes();
//            TelePrompter::printTree(1, "Gradient", nNodes, time);
//        }
//        timer.resume();
//        T *result = new T;
//        result->add(trees, 0);
//        double time = timer.elapsed();
//        int nNodes = result->getNNodes();
//        TelePrompter::printTree(1, "Sum divergence", nNodes, time);
//        delete trees[0];
//        delete trees[1];
//        delete trees[2];
//        return result;
//    }

    template<class T>
    T** allocPtrArray(int n_funcs) {
        T **ptr = new T*[n_funcs];
        for (int i = 0; i < n_funcs; i++) {
            ptr[i] = 0;
        }
        return ptr;
    }

    template<class T>
    void clearPtrArray(int n_funcs, T **ptr) {
        if (ptr == 0) MSG_FATAL("Clearing NULL pointer");
        for (int i = 0; i < n_funcs; i++) {
            if (ptr[i] != 0) {
                delete ptr[i];
            }
            ptr[i] = 0;
        }
    }

    template<class T>
    T** deletePtrArray(int n_funcs, T ***ptr) {
        if (*ptr != 0) {
            for (int i = 0; i < n_funcs; i++) {
                if ((*ptr)[i] != 0) {
                    delete (*ptr)[i];
                }
                (*ptr)[i] = 0;
            }
            delete[] *ptr;
        }
        return 0;
    }
};

#endif // XCOPERATOR_H

