#pragma once

#include "MRCPP/Parallel"
#include <map>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

namespace omp {
extern int n_threads;
} // namespace omp

namespace mpi {
extern bool numerically_exact;
extern int shared_memory_size;

extern int world_rank;
extern int world_size;
extern int orb_rank;
extern int orb_size;
extern int share_rank;
extern int share_size;
extern int sh_group_rank;
extern int is_bank;
extern int is_bankclient;
extern int bank_size;
extern std::vector<int> bankmaster;

extern MPI_Comm comm_orb;
extern MPI_Comm comm_share;
extern MPI_Comm comm_sh_group;
extern MPI_Comm comm_bank;

void initialize();
void finalize();
void barrier(MPI_Comm comm);

bool grand_master();
bool share_master();
bool my_orb(const Orbital &orb);
bool my_unique_orb(const Orbital &orb);
void distribute(OrbitalVector &Phi);
void free_foreign(OrbitalVector &Phi);
OrbitalChunk get_my_chunk(OrbitalVector &Phi);

void send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm = mpi::comm_orb);
void recv_orbital(Orbital &orb, int src, int tag, MPI_Comm comm = mpi::comm_orb);

void send_function(QMFunction &func, int dst, int tag, MPI_Comm comm);
void recv_function(QMFunction &func, int src, int tag, MPI_Comm comm);
void share_function(QMFunction &func, int src, int tag, MPI_Comm comm);

void reduce_function(double prec, QMFunction &func, MPI_Comm comm);
void broadcast_function(QMFunction &func, MPI_Comm comm);

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

struct deposit {
    Orbital *orb;
    int id;     // to identify what is deposited
    int source; // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};

class Bank {
public:
    Bank();
    ~Bank();
    int open();
    int close();
    int clear_all(int i, MPI_Comm comm);
    int clear(int ix);
    int put_orb(int id, Orbital &orb);
    int get_orb(int id, Orbital &orb);
    int put_func(int id, QMFunction &func);
    int get_func(int id, QMFunction &func);

private:
    int const CLOSE_BANK = 1;
    int const CLEAR_BANK = 2;
    int const GET_ORBITAL = 3;
    int const SAVE_ORBITAL = 4;
    int const GET_FUNCTION = 5;
    int const SAVE_FUNCTION = 6;
    int clear_bank();
    std::map<int, int> id2ix;
    std::vector<deposit> deposits;
    std::map<int, int> id2qu;
    std::vector<queue_struct> queue;
};

namespace mpi {
extern Bank orb_bank;
} // namespace mpi

} // namespace mrchem
