#include "MREnv.h"
#include "Timer.h"
#include "mrchem.h"

using namespace std;

void MREnv::initializeMRCPP(int argc, char **argv) {
    int nThreads = omp_get_max_threads();
    int nHosts = 1;

    omp_set_dynamic(0);
    Eigen::setNbThreads(1);

    const char *infile = 0;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    Input = Getkw(infile, false, true);

    int printlevel = Input.get<int>("printlevel");
    bool teletype = Input.get<bool>("teletype");

    TelePrompter::init(printlevel, teletype, "MRCHEM");
    TelePrompter::setPrecision(15);

    println(0,endl << endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***    MRChem " << PROGRAM_VERSION << " (rev. " <<
            GIT_REVISION << ")                       ***");
    println(0,"***                                                      ***");
    println(0,"***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0,"***    Jonas Juselius   <jonas.juselius@uit.no>          ***");
    println(0,"***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);
    println(0,"Print level  : " <<  printlevel << endl);

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

void MREnv::finalizeMRCPP(Timer t) {
    double wt = t.getWallTime();
    double ut = t.getUserTime();
    double st = t.getSystemTime();
    SET_PRINT_PRECISION(6);
    println(0,endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***                    Exiting MRChem                    ***");
    println(0,"***                                                      ***");
    println(0,"***                 Wall:   " << wt << "                 ***");
    println(0,"***                 User:   " << ut << "                 ***");
    println(0,"***                 System: " << st << "                 ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);
}
