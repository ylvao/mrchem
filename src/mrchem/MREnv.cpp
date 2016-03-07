#include "MREnv.h"

using namespace std;

void MREnv::initializeMRCPP() {
    int nThreads = omp_get_max_threads();
    int nHosts = node_group.size();

    omp_set_dynamic(0);
    Eigen::setNbThreads(1);

    int printLevel = 10;
    bool teletype = true;

    TelePrompter::init(printLevel, teletype, "MRChem");

    println(0,endl << endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***    MRChem " << PROJECT_VERSION << " (rev. " <<
            GIT_REVISION << ")                       ***");
    println(0,"***                                                      ***");
    println(0,"***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0,"***    Jonas Juselius   <jonas.juselius@uit.no>          ***");
    println(0,"***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);
    println(0,"Print level  : " <<  printLevel << endl);

    if (nHosts > 1 or nThreads > 1) {
        println(0,"+++ Parallel execution: ");
        println(0,"  MPI hosts available     : " << nHosts);
        println(0,"  Threads/host            : " << nThreads);
        println(0,"  Total used CPUs         : " << nHosts*nThreads);
        println(0,"");
    } else {
        println(0,"+++ Serial execution" << endl);
    }
}

void MREnv::finalizeMRCPP(double t) {
    SET_PRINT_PRECISION(5);
    println(0,endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***                     Exiting MRCPP                    ***");
    println(0,"***                                                      ***");
    println(0,"***               World clock: " << t << "               ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);
}
