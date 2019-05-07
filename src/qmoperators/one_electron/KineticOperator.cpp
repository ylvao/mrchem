#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "KineticOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Expectation value matrix
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix KineticOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    mrcpp::print::header(1, "Compute Kinetic Matrix Elements");

    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T_x = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_y = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_z = ComplexMatrix::Zero(Ni, Nj);
    {
        Timer timer;
        if (&bra == &ket) {
            OrbitalVector dKet = p_x(ket);
            T_x = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_x(bra);
            OrbitalVector dKet = p_x(ket);
            T_x = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::time(1, "T_x", timer);
    }
    {
        Timer timer;
        if (&bra == &ket) {
            OrbitalVector dKet = p_y(ket);
            T_y = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_y(bra);
            OrbitalVector dKet = p_y(ket);
            T_y = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::time(1, "T_y", timer);
    }
    {
        Timer timer;
        if (&bra == &ket) {
            OrbitalVector dKet = p_z(ket);
            T_z = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_z(bra);
            OrbitalVector dKet = p_z(ket);
            T_z = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::time(1, "T_z", timer);
    }
    mrcpp::print::footer(1, timer, 2);

    return 0.5 * (T_x + T_y + T_z);
}

/** @brief Expectation value (dagger version)
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * NOT IMPLEMENTED
 */
ComplexMatrix KineticOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
