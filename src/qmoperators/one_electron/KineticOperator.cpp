#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "KineticOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/Bank.h"

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
    mrcpp::print::memory(2, "memusage Kinetic start");

    // we keep both methods for testing purposes
    if (mpi::bank_size > 0) {
        Timer timer;
        // 1) save all derivatives
        OrbitalVector dKet = p_x(ket);
        mrcpp::print::tree(2, "compute xderivatives ", 1, 1, timer.elapsed());
        mpi::barrier(mpi::comm_orb);
        for (int i = 0; i < bra.size(); i++) {
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
        }
        mrcpp::print::tree(2, "save xderivatives ", 1, 1, timer.elapsed());
        // 2) symmetric dot product
        mpi::barrier(mpi::comm_orb);
        mrcpp::print::tree(2, "barrier ", 1, 1, timer.elapsed());
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
            //            mrcpp::print::tree(2, "get orbi ", i, 1, timer.elapsed());
            for (int j = 0; j < ket.size(); j++) {
                if (not mpi::my_orb(ket[j])) continue;
                // if( (&bra == &ket) and (i>(Ni-1)/2+j or (i>=j-(Ni-1)/2 and i<j))) continue;
                if ((&bra == &ket) and (i + j + (i > j)) % 2) continue;
                Orbital derj;
                if (&bra == &ket) {
                    derj = dKet[j];
                } else {
                    derj = p_z(ket[j]);
                }
                //                mrcpp::print::tree(2, "get bankj ", i, j, timer.elapsed());

                T_x(i, j) = orbital::dot(deri, derj);
                // mrcpp::print::tree(2, "get dot ", i, j, timer.elapsed());
                if (&bra == &ket) T_x(j, i) = std::conj(T_x(i, j));
            }
            mrcpp::print::tree(2, "j products ", i, 1, timer.elapsed());
        }
        mrcpp::print::tree(2, "Tx ready ", 1, 1, timer.elapsed());
        mpi::barrier(mpi::comm_orb);
        mrcpp::print::tree(2, "barrier ", 1, 1, timer.elapsed());

        dKet = p_y(ket);
        mpi::barrier(mpi::comm_orb);
        for (int ii = shift; ii < bra.size() + shift; ii++) {
            int i = ii % bra.size();
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
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
                //                if( (&bra == &ket) and (i>(Ni-1)/2+j or (i>=j-(Ni-1)/2 and i<j))) continue;
                if ((&bra == &ket) and (i + j + (i > j)) % 2) continue;
                Orbital derj;
                if (&bra == &ket) {
                    derj = dKet[j];
                } else {
                    derj = p_z(ket[j]);
                }
                //       mrcpp::print::tree(2, "get bankj ", i, j, timer.elapsed());
                T_y(i, j) = orbital::dot(deri, derj);
                // mrcpp::print::tree(2, "get dot ", i, j, timer.elapsed());
                if (&bra == &ket) T_y(j, i) = std::conj(T_y(i, j));
            }
        }
        dKet = p_z(ket);
        mpi::barrier(mpi::comm_orb);
        for (int i = 0; i < bra.size(); i++) {
            if (not mpi::my_orb(bra[i])) continue;
            Orbital der = dKet[i];
            mpi::orb_bank.put_orb(i, der);
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
        T_tot = 0.5 * (T_x + T_y + T_z);
        mpi::allreduce_matrix(T_tot, mpi::comm_orb);
        mrcpp::print::memory(2, "memusage Kinetic end1");
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
        mrcpp::print::memory(2, "memusage Kinetic end2");
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
