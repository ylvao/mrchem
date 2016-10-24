#include "MREnv.h"
#include "Timer.h"
#include "mrchem.h"
#include "LegendreBasis.h"
#include "InterpolatingBasis.h"

using namespace std;

MultiResolutionAnalysis<3> *MRA;

void MREnv::initializeMRCPP(int argc, char **argv) {
#ifdef HAVE_MPI
    MPI_Init(NULL, NULL);
    MPI_Initializations();
#endif
    int nThreads = omp_get_max_threads();

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

    if (MPI_size > 1 or nThreads > 1) {
        println(0,"+++ Parallel execution: ");
        println(0,"  MPI hosts available     : " << MPI_size);
        println(0,"  Threads/host            : " << nThreads);
        println(0,"  Total used CPUs         : " << MPI_size*nThreads);
        println(0,"");
    } else {
        println(0,"+++ Serial execution" << endl);
    }

    // Initialize global MRA
    initializeMRA();
}

void MREnv::finalizeMRCPP(const Timer t) {
    // Delete global MRA
    if (MRA != 0) delete MRA;
    MRA = 0;

    double wt = t.getWallTime();
    SET_PRINT_PRECISION(6);
    println(0,endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***                    Exiting MRChem                    ***");
    println(0,"***                                                      ***");
    println(0,"***               Wall time: " << wt << "                ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

void MREnv::initializeMRA() {
    // Constructing world box
    int scale = Input.get<int>("World.scale");
    vector<int> corner = Input.getIntVec("World.corner");
    vector<int> boxes = Input.getIntVec("World.boxes");
    NodeIndex<3> idx(scale, corner.data());
    BoundingBox<3> world(idx, boxes.data());

    // Constructing scaling basis
    int order = Input.get<int>("order");
    string wtype = Input.get<string>("wavelet");

    // Initializing MRA
    if (wtype == "I") {
        InterpolatingBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(world, basis);
    } else if (wtype == "L") {
        LegendreBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(world, basis);
    } else {
        MSG_FATAL("Invalid wavelet type!");
    }
}
