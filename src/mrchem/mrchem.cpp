/** \mainpage The MRChem program
 *
 * \author Stig Rune Jensen
 *
 * \version 1.0
 *
 * \par Copyright:
 * GPLv4
 *
 */

#include <Eigen/Core>

#include "MREnv.h"

int main(int argc, char **argv) {
    mpi::environment env(argc, argv);

    Timer rolex;
    rolex.restart();

    MREnv::initializeMRCPP(argc, argv);
    MREnv::finalizeMRCPP(rolex);

    return 0;
}

