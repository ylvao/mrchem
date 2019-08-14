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
    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T_x = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_y = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_z = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_tot = ComplexMatrix::Zero(Ni, Nj);

    // we keep both methods for testing purposes
    if (mpi::bank_size > 0) {
        Timer timer;
        // 1) save all derivatives
        OrbitalVector dKet = p_x(ket);
        int nNodes = 0, sNodes = 0;
        mpi::barrier(mpi::comm_orb);
        for (int i = 0; i < bra.size(); i++) {
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
            nNodes += der.getNNodes(NUMBER::Total);
            sNodes += der.getSizeNodes(NUMBER::Total);
        }
        // 2) symmetric dot product
        mpi::barrier(mpi::comm_orb);
        int shift = mpi::orb_rank;
        for (int ii = shift; ii < bra.size() + shift; ii++) {
            int i = ii % bra.size();
            Orbital deri;
            if (mpi::my_orb(bra[i])) {
                deri = dKet[i];
            } else {
                Timer timer;
                mpi::orb_bank.get_orb(i, deri);
            }
            for (int j = 0; j < ket.size(); j++) {
                if (not mpi::my_orb(ket[j])) continue;
                if ((&bra == &ket) and (i + j + (i > j)) % 2) continue;
                Orbital derj;
                if (&bra == &ket) {
                    derj = dKet[j];
                } else {
                    derj = p_z(ket[j]);
                }
                T_x(i, j) = orbital::dot(deri, derj);
                if (&bra == &ket) T_x(j, i) = std::conj(T_x(i, j));
            }
        }
        mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
        timer.start();
        mpi::barrier(mpi::comm_orb);
        nNodes = 0, sNodes = 0;
        dKet = p_y(ket);
        mpi::barrier(mpi::comm_orb);
        for (int ii = shift; ii < bra.size() + shift; ii++) {
            int i = ii % bra.size();
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
            nNodes += der.getNNodes(NUMBER::Total);
            sNodes += der.getSizeNodes(NUMBER::Total);
        }
        mpi::barrier(mpi::comm_orb);
        // 2) symmetric dot product
        for (int i = 0; i < bra.size(); i++) {
            Orbital deri;
            if (mpi::my_orb(bra[i])) {
                deri = dKet[i];
            } else {
                mpi::orb_bank.get_orb(i, deri);
            }
            for (int j = 0; j < ket.size(); j++) {
                if (not mpi::my_orb(ket[j])) continue;
                if ((&bra == &ket) and (i + j + (i > j)) % 2) continue;
                Orbital derj;
                if (&bra == &ket) {
                    derj = dKet[j];
                } else {
                    derj = p_z(ket[j]);
                }
                T_y(i, j) = orbital::dot(deri, derj);
                if (&bra == &ket) T_y(j, i) = std::conj(T_y(i, j));
            }
        }
        mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
        timer.start();
        dKet = p_z(ket);
        mpi::barrier(mpi::comm_orb);
        nNodes = 0, sNodes = 0;
        for (int i = 0; i < bra.size(); i++) {
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
            nNodes += der.getNNodes(NUMBER::Total);
            sNodes += der.getSizeNodes(NUMBER::Total);
        }
        mpi::barrier(mpi::comm_orb);
        // 2) symmetric dot product
        for (int ii = shift; ii < bra.size() + shift; ii++) {
            int i = ii % bra.size();
            Orbital deri;
            if (mpi::my_orb(bra[i])) {
                deri = dKet[i];
            } else {
                mpi::orb_bank.get_orb(i, deri);
            }
            for (int j = 0; j < ket.size(); j++) {
                if (not mpi::my_orb(ket[j])) continue;
                if ((&bra == &ket) and (i + j + (i > j)) % 2) continue;
                Orbital derj;
                if (&bra == &ket) {
                    derj = dKet[j];
                } else {
                    derj = p_z(ket[j]);
                }
                T_z(i, j) = orbital::dot(deri, derj);
                if (&bra == &ket) T_z(j, i) = std::conj(T_z(i, j));
            }
        }
        mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
        T_tot = 0.5 * (T_x + T_y + T_z);
        mpi::allreduce_matrix(T_tot, mpi::comm_orb);
        mpi::orb_bank.clear_all(mpi::orb_rank, mpi::comm_orb);
        return T_tot;
    } else {
        {
            Timer timer;
            int nNodes = 0, sNodes = 0;
            if (&bra == &ket) {
                OrbitalVector dKet = p_x(ket);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dKet);
                T_x = orbital::calc_overlap_matrix(dKet);
            } else {
                OrbitalVector dBra = p_x(bra);
                OrbitalVector dKet = p_x(ket);
                nNodes += orbital::get_n_nodes(dBra);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dBra);
                sNodes += orbital::get_size_nodes(dKet);
                T_x = orbital::calc_overlap_matrix(dBra, dKet);
            }
            mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
        }
        {
            Timer timer;
            int nNodes = 0, sNodes = 0;
            if (&bra == &ket) {
                OrbitalVector dKet = p_y(ket);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dKet);
                T_y = orbital::calc_overlap_matrix(dKet);
            } else {
                OrbitalVector dBra = p_y(bra);
                OrbitalVector dKet = p_y(ket);
                nNodes += orbital::get_n_nodes(dBra);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dBra);
                sNodes += orbital::get_size_nodes(dKet);
                T_y = orbital::calc_overlap_matrix(dBra, dKet);
            }
            mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
        }
        {
            Timer timer;
            int nNodes = 0, sNodes = 0;
            if (&bra == &ket) {
                OrbitalVector dKet = p_z(ket);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dKet);
                T_z = orbital::calc_overlap_matrix(dKet);
            } else {
                OrbitalVector dBra = p_z(bra);
                OrbitalVector dKet = p_z(ket);
                nNodes += orbital::get_n_nodes(dBra);
                nNodes += orbital::get_n_nodes(dKet);
                sNodes += orbital::get_size_nodes(dBra);
                sNodes += orbital::get_size_nodes(dKet);
                T_z = orbital::calc_overlap_matrix(dBra, dKet);
            }
            mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
        }
        return 0.5 * (T_x + T_y + T_z);
    }
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
