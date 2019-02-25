/** \mainpage The MRChem main program
 *
 * \author Stig Rune Jensen
 *
 * \version 1.0
 *
 * \par Copyright:
 * GPLv4
 *
 */

#include "MRCPP/Timer"

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

#include "SCFDriver.h"

nlohmann::json mrchem::json_input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    mpi::initialize(argc, argv);
    mrenv::initialize(argc, argv);

    Timer timer;

    /*
    SCFDriver driver(input);
    driver.setup();
    driver.run();
    driver.clear();
    */

    timer.stop();
    double wt = timer.getWallTime();

    mrenv::finalize(wt);
    mpi::finalize();
    return 0;
}
