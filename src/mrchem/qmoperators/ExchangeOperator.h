#ifndef EXCHANGEOPERATOR_H
#define EXCHANGEOPERATOR_H

#include <Eigen/Core>

#include "QMOperator.h"
#include "OrbitalAdder.h"
#include "OrbitalMultiplier.h"
#include "GridCleaner.h"
#include "GridGenerator.h"
#include "PoissonOperator.h"

class OrbitalVector;
template<int D> class MultiResolutionAnalysis;

class ExchangeOperator : public QMOperator {
public:
    ExchangeOperator(double build_prec,
                     const MultiResolutionAnalysis<3> &mr,
                     OrbitalVector &phi,
                     double x_fac);
    virtual ~ExchangeOperator();

    virtual void setup(double prec) = 0;
    virtual void clear() = 0;

    void setExchangeFactor(double x_fac) { this->x_factor = x_fac; }
    double getExchangeFactor() const { return this->x_factor; }

    void setScreen(bool s) { this->screen = s; }
    bool getScreen() const { return this->screen; }

protected:
    OrbitalAdder add;
    OrbitalMultiplier mult;
    GridGenerator<3> grid;
    PoissonOperator poisson;   ///< Poisson operator to compute potential

    double x_factor;            ///< Exchange factor for Hybrid XC functionals
    OrbitalVector *orbitals_0;  ///< The orbitals that define the exchange

    bool screen;                ///< Apply screening in exchange evaluation
    Eigen::VectorXd tot_norms;  ///< Total norms for use in screening
    Eigen::MatrixXd part_norms; ///< Partial norms for use in screening

    double getScaledPrecision(int i, int j) const;
};

#endif // EXCHANGEOPERATOR_H
