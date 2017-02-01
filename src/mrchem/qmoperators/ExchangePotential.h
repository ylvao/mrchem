#ifndef EXCHANGEPOTENTIAL_H
#define EXCHANGEPOTENTIAL_H

#include "ExchangeOperator.h"
#include "OrbitalVector.h"

/*! \file ExchangePotential.h
 *  \class ExchangePotential
 *  \brief Hartree-Fock exchange potential for a given orbital.
 *  \author Stig Rune Jensen
 *  \date 2015
 *
 */

class ExchangePotential : public ExchangeOperator {
public:
    
  /*! Constructor for the Exchange Potential 
   * \param[in] prec prcision requested to the calculation
   * \param[in] phi vector of orbitals which define the exchange operator
   * \param[in] x_fac Exchange factor for Hybrid XC functionals
   */
    ExchangePotential(double prec, OrbitalVector &phi, double x_fac = 1.0);
    virtual ~ExchangePotential() { }

    /*! Rotates the Exchange Fock Matrix
     * \param[in] U unitary matrix defining the rotation
     */
    virtual void rotate(Eigen::MatrixXd &U);

    /*! Set the precision for the application of the Exchange operator
     * \param[in] prec reqested precision
     */
    virtual void setup(double prec);

    /*! Clears the Exchange Operator
     *  Clears both the set of orbitals used to define the Exchange Potential and the underlying Exchange operator.
     */
    virtual void clear();

    /*! Applies operator potential
     *  \param[in] phi_p input orbital
     *
     * The exchange potential is applied to the given orbital. 
     */
    virtual Orbital* operator() (Orbital &phi_p);

    /*! Applies the adjoint of the operator
     *  \param[in] phi_p input orbital
     *
     * NOT IMPLEMENTED
     */
    virtual Orbital* adjoint(Orbital &phi_p);

    using QMOperator::operator();
    using QMOperator::adjoint;
protected:
    OrbitalVector exchange_0;  ///< Set to keep precomputed exchange orbitals

    /*! computes the exchange potential
     *  \param[in] phi_p input orbital
     *
     * The exchange potential is computed and applied on the fly to the given orbital. 
     */
    Orbital* calcExchange(Orbital &phi_p);
    /*! Test if a given contribution has been precomputed
     * \param[in] phi_p orbital for which the check is performed
     *
     * If the given contribution has been precomputed, it is simply added, without additional recalculation
     */
    Orbital* testPreComputed(Orbital &phi_p);

    /*! precomputes the exchange potential
     *  \param[in] phi_p input orbital
     *
     * The exchange potential is (pre)computed but not applied to any orbital. 
     */
    void calcInternalExchange();

    /*! computes the diagonal part of the exchange potential
     *  \param[in] i orbital index
     *
     * The diagonal term X_ii is computed. 
     */
    int calcInternal(int i);
    /*! computes the off-diagonal part of the exchange potential
     *  \param[in] i first orbital index
     *  \param[in] j second orbital index
     *
     * The off-diagonal term X_ij is computed. 
     */
    int calcInternal(int i, int j);
};

#endif // EXCHANGEPOTENTIAL_H
