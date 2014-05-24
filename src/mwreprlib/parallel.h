/*
 *
 *
 *  \date Sep 27, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Convenience definitions for parallel builds
 */

#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <boost/serialization/serialization.hpp>

#include "config.h"

int get_locale_index_range(int rank, int nWork, int &start, int &end);
int get_locale_index(int nWork, int idx);
bool locale_needs_sync(int nWork);
int get_hosts_bitmap_size();

#ifdef HAVE_OPENMP

#define OPENMP
#include <omp.h>

#else

#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_set_dynamic(n)
#define omp_set_lock(x)
#define omp_unset_lock(x)
#define omp_test_lock(x)

#endif

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#include <boost/mpi/timer.hpp>
namespace mpi = boost::mpi;
#define BOOST_MPI_HOMOGENEOUS

/** Return the number MPI hosts available. */
inline int get_mpi_world_size() {
	mpi::communicator world;
	return world.size();
}

#else
namespace mpi {
	struct communicator {
		int rank() const { return 0; }
		int size() const { return 1; }
		void barrier() const { }
	};
	struct environment {
		environment(int argc, char **argv) { }
	};
	struct timer {
		void restart() { }
		double elapsed() { return 0.0; }
	};
	typedef int request;
}
inline int get_mpi_world_size() { return 0; }
#endif

#endif /* PARALLEL_H_ */
