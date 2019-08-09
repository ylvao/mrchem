/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "utils/Bank.h"
#include "MRCPP/Timer"

// TODO: bank for serial runs

int const func_id_shift = 10000; // to ensure that orbitals and function do not collide
using mrcpp::Timer;

namespace mrchem {

Bank::Bank() {}
Bank::~Bank() {
    for (int ix = 1; ix < this->deposits.size(); ix++) this->clear(ix);
}

int Bank::open() {
    std::cout << " memory bank opened " << std::endl;
    bool listen = true;
    MPI_Status status;
    char safe_data1;
    int deposit_size = sizeof(deposit);
    int n_chunks, ix;
    int message;
    // TODO: test queue for incoming requests

    deposits.resize(1); // we reserve 0, since it is the value returned for undefined key
    queue.resize(1);    // we reserve 0, since it is the value returned for undefined key

    // The bank never goes out of this loop as long as it receives a close signal!
    while (listen) {
        MPI_Recv(&message, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, mpi::comm_bank, &status);
        //        std::cout << "world_rank " << mpi::world_rank << " message " << message << " Bank got request from
        //        client "
        //          << status.MPI_SOURCE << " id " << status.MPI_TAG << std::endl;
        if (message == CLOSE_BANK) {
            std::cout << "Bank is closing" << std::endl;
            break; // close bank. To be made cleaner
        }
        if (message == CLEAR_BANK) {
            std::cout << "Bank removes all deposits" << std::endl;
            this->clear_bank();
            // send message that it is ready (value of message is not used)
            MPI_Send(&message, 1, MPI_INTEGER, status.MPI_SOURCE, 77, mpi::comm_bank);
            break; //
        }
        if (message == GET_ORBITAL or message == GET_FUNCTION) {
            // withdrawal
            int ix = id2ix[status.MPI_TAG];
            if (ix == 0) {
                // the id does not exist. Put in queue and Wait until it is defined
                std::cout << ix << " requext id " << status.MPI_TAG << " from " << status.MPI_SOURCE << " is QUEUED"
                          << std::endl;
                if (id2qu[status.MPI_TAG] == 0) {
                    queue.push_back({status.MPI_TAG, {status.MPI_SOURCE}});
                    id2qu[status.MPI_TAG] = queue.size() - 1;
                } else {
                    // somebody is already waiting for this id. queue in queue
                    queue[id2qu[status.MPI_TAG]].clients.push_back(status.MPI_SOURCE);
                }
            } else {
                if (deposits[ix].id != status.MPI_TAG) std::cout << ix << " Bank acounting error " << std::endl;
                //                std::cout << ix << " memory bank will send orbital " << status.MPI_TAG << " to " <<
                //                status.MPI_SOURCE
                //        << std::endl;
                if (message == GET_ORBITAL) {
                    Timer timer;
                    mpi::send_orbital(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, mpi::comm_bank);
                    // std::cout<<deposits[ix].id<<" sent from bank in "<< (int)(1000*timer.elapsed())<<"
                    // "<<deposits[ix].orb->getNNodes(0)<<std::endl;
                }
                if (message == GET_FUNCTION) {
                    mpi::send_function(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, mpi::comm_bank);
                }
            }
        }
        if (message == SAVE_ORBITAL or message == SAVE_FUNCTION) {
            // make a new deposit
            if (id2ix[status.MPI_TAG]) {
                //  std::cout << "WARNING: " << deposits[ix].id << " already deposited; will be overwritten" <<
                //  std::endl;
                ix = id2ix[status.MPI_TAG];
                //                this->clear(ix);
            } else {
                ix = deposits.size(); // NB: ix is now index of last element + 1
                deposits.resize(ix + 1);
            }
            deposits[ix].orb = new Orbital(0);
            // std::cout << "world_rank " << mpi::world_rank << " " << ix << " Bank new deposit id " << status.MPI_TAG
            //        << " from " << status.MPI_SOURCE << std::endl;

            deposits[ix].id = status.MPI_TAG;
            id2ix[deposits[ix].id] = ix;
            deposits[ix].source = status.MPI_SOURCE;
            if (message == SAVE_ORBITAL) {
                mpi::recv_orbital(*deposits[ix].orb, deposits[ix].source, status.MPI_TAG, mpi::comm_bank);
            }
            if (message == SAVE_FUNCTION) {
                mpi::recv_function(*deposits[ix].orb, deposits[ix].source, status.MPI_TAG, mpi::comm_bank);
            }
            //            std::cout << mpi::world_rank << " " << ix << " Bank got id " << status.MPI_TAG << " from "
            //          << status.MPI_SOURCE << std::endl;
            if (id2qu[deposits[ix].id] != 0) {
                // someone is waiting for those data. Send to them
                int iq = id2qu[deposits[ix].id];
                if (deposits[ix].id != queue[iq].id) std::cout << ix << " Bank queue acounting error " << std::endl;
                for (int iqq : queue[iq].clients) {
                    std::cout << ix << " memory bank will send orbital (queued) " << queue[iq].id << " to " << iqq
                              << std::endl;
                    if (message == SAVE_ORBITAL) {
                        mpi::send_orbital(*deposits[ix].orb, iqq, queue[iq].id, mpi::comm_bank);
                    }
                    if (message == SAVE_FUNCTION) {
                        mpi::send_function(*deposits[ix].orb, iqq, queue[iq].id, mpi::comm_bank);
                    }
                }
                //                for (auto &a : queue) { std::cout << a.id << " queue before erase " << std::endl; }
                queue.erase(queue.begin() + iq);
                // for (auto &a : queue) { std::cout << a.id << " queue after erase " << std::endl; }
                id2qu.erase(deposits[ix].id);
            }
        }
    }
}

// save orbital in Bank with identity id
int Bank::put_orb(int id, Orbital &orb) {
    // for now we distribute according to id
    MPI_Send(&SAVE_ORBITAL, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    //    std::cout << mpi::orb_rank << " deposited id " << id << std::endl;
    mpi::send_orbital(orb, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
}

// get orbital with identity id
int Bank::get_orb(int id, Orbital &orb) {
    MPI_Status status;

    MPI_Send(&GET_ORBITAL, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    Timer t;
    mpi::recv_orbital(orb, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    // std::cout << mpi::orb_rank << " received " << id <<" after "<<std::setprecision(2) << std::fixed
    // <<1000*t.elapsed() <<" msec"<< std::endl;
}

// save function in Bank with identity id
int Bank::put_func(int id, QMFunction &func) {
    // for now we distribute according to id
    id += func_id_shift;
    MPI_Send(&SAVE_FUNCTION, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    //    std::cout << mpi::orb_rank << " deposited id " << id << std::endl;
    mpi::send_function(func, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
}

// get function with identity id
int Bank::get_func(int id, QMFunction &func) {
    MPI_Status status;
    MPI_Send(&GET_FUNCTION, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    id += func_id_shift;
    mpi::recv_function(func, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    //    std::cout << mpi::orb_rank << " received " << id << std::endl;
}

// Ask to close the Bank
int Bank::close() {
    for (int i = 0; i < mpi::bank_size; i++) {
        MPI_Send(&CLOSE_BANK, 1, MPI_INTEGER, mpi::bankmaster[i], 0, mpi::comm_bank);
    }
}

// remove all deposits
// NB:: collective call. All clients must call this
int Bank::clear_all(int i, MPI_Comm comm) {
    // 1) wait until all clients are ready
    mpi::barrier(comm);
    // master send signal to bank
    if (i == 0) {
        for (int i = 0; i < mpi::bank_size; i++) {
            MPI_Send(&CLEAR_BANK, 1, MPI_INTEGER, mpi::bankmaster[i], 0, mpi::comm_bank);
            ;
        }
        for (int i = 0; i < mpi::bank_size; i++) {
            // wait until Bank is finished and has sent signal
            MPI_Status status;
            int message;
            MPI_Recv(&message, 1, MPI_INTEGER, mpi::bankmaster[i], 77, mpi::comm_bank, &status);
        }
    }
    mpi::barrier(comm);
}

int Bank::clear_bank() {
    for (int ix = 1; ix < this->deposits.size(); ix++) this->clear(ix);
    this->deposits.resize(1);
    this->queue.resize(1);
    this->id2ix.clear();
    this->id2qu.clear();
}

int Bank::clear(int ix) {
    delete deposits[ix].orb;
}

} // namespace mrchem
