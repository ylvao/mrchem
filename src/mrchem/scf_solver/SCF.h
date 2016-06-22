#ifndef SCF_H
#define SCF_H

#include <vector>
#include <Eigen/Core>

#include "OrbitalAdder.h"
#include "TelePrompter.h"

class HelmholtzOperatorSet;
class Accelerator;
class FockOperator;
class OrbitalVector;
class Orbital;

class SCF {
public:
    SCF(const MultiResolutionAnalysis<3> &mra, HelmholtzOperatorSet &h);
    virtual ~SCF();

    double getOrbitalPrecision() const { return this->orbPrec[0]; }
    double getOrbitalThreshold() const { return this->orbThrs; }
    double getPropertyThreshold() const { return this->propThrs; }

    void setThreshold(double orb_thrs, double prop_thrs);
    void setOrbitalPrec(double init, double final);
    void setMaxIterations(int m_iter) { this->maxIter = m_iter; }

    void setRotation(int iter) { this->iterPerRotation = iter; }
    void setOrthogonalize() { this->iterPerRotation = 0; }
    void setDiagonalize(int iter) { this->iterPerRotation = -abs(iter); }
    void setLocalize(int iter) { this->iterPerRotation = abs(iter); }

protected:
    int nIter;
    int maxIter;
    int iterPerRotation; ///< Positive localization, negative diagonalization, zero S_m12
    double orbThrs;  ///< Convergence threshold orbital update norm
    double propThrs; ///< Convergence threshold property
    double orbPrec[3];

    std::vector<double> orbError;
    std::vector<double> property;

    HelmholtzOperatorSet *helmholtz;// Pointer to external object, do not delete!
    OrbitalAdder add;
//    OrbitalRotator rotate;

    bool needLocalization() const;
    bool needDiagonalization() const;

    void adjustPrecision(double error);
    void resetPrecision();

    void printUpdate(const std::string &name, double P, double dP) const;
    double getUpdate(const std::vector<double> &vec, int i, bool absPrec) const;

    void printOrbitals(const Eigen::MatrixXd &F, const OrbitalVector &phi) const;
    void printConvergence(bool converged) const;
    void printCycle() const;
    void printTimer(double t) const;
    void printMatrix(int level, const Eigen::MatrixXd &M,
                     const char &name, int pr = 5) const;

    bool accelerate(Accelerator *acc,
                    OrbitalVector *phi,
                    OrbitalVector *d_phi,
                    Eigen::MatrixXd *F = 0,
                    Eigen::MatrixXd *dF = 0);

    void applyHelmholtzOperators(OrbitalVector &phi_np1,
                                 OrbitalVector &phi_n,
                                 Eigen::MatrixXd &F_n,
                                 bool adjoint = false);

    virtual Orbital* getHelmholtzArgument(int i,
                                          OrbitalVector &phi,
                                          Eigen::MatrixXd &F,
                                          bool adjoint) = 0;

    Orbital* calcMatrixPart(int i,
                            Eigen::MatrixXd &M,
                            OrbitalVector &phi);
};

#endif // SCF_H

