#include "MREnv.h"
#include "TelePrompter.h"

#include "QuadratureCache.h"

using namespace std;

void MREnv::initializeMRCPP(int k, double prec) {
    omp_set_dynamic(0);
    mpi::communicator world;

    int nThreads = omp_get_max_threads();
    Eigen::setNbThreads(1); 

    int printLevel = 0;
    bool teletype = 0;

    TelePrompter::init(printLevel, teletype, "MRGRID");

    println(0, endl << endl);
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***    MRGrid " << PROJECT_VERSION << " (rev. " <<
	                        GIT_REVISION << ")                       ***");
    println(0, "***                                                      ***");
    println(0, "***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0, "***    Jonas Juselius   <jonas.juselius@uit.no>          ***");
    println(0, "***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, endl);
    println(0, "Print level  : " <<  printLevel << endl);

    if (world.size() > 1 or nThreads > 1) {
	println(0, "+++ Parallel execution: ");
	println(0, "  Num MPI hosts available : " << world.size());
	println(0, "  Threads/host            : " << nThreads);
	println(0, "  Total used CPUs         : " << world.size()*nThreads);
	println(0, "");
    } else {
	println(0, "+++ Serial execution" << endl);
    }

    // initialize QuadratureCache globally to [0.1]
    getQuadratureCache(qCache);
    qCache.setBounds(0.0, 1.0);

    //Initialize world
    int uni_depth = 0;
    int max_depth = 20;
    int polytype = Interpol;

    int root_scale    = -5;
    int nbox[3]      = {    1,    1,    1};
    int transl[3]    = {    0,    0,    0};
    double origin[3] = {-32.0,-32.0,-32.0};

    //BoundingBox<1>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    //BoundingBox<2>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    //BoundingBox<3>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());

    //const BoundingBox<3> &worldbox = BoundingBox<3>::getWorldBox();
    //println(0, worldbox);

    println(0, "+++ Parameters:");
    println(0, "  Initial scale:    " << root_scale);
    println(0, "  Uniform depth:    " << uni_depth);
    println(0, "  Max depth:        " << max_depth);
    println(0, "  Precision:        " << prec);
    println(0, "  Polynomial order: " << k);
    printout(0, "  Polynomial type:  ");
    if (polytype == Interpol) println(0, "Interpolating");
    if (polytype == Legendre) println(0, "Legendre");
    println(0, endl);

    //initializeTrees(order, max_depth, rel_prec, polytype);
}

void MREnv::initializeTrees(int k, int depth, double prec, int type) {
/*
    FunctionTree<1>::setDefaultOrder(k);
    FunctionTree<1>::setDefaultMaxDepth(depth);
    FunctionTree<1>::setDefaultPrecision(prec);
    FunctionTree<1>::setDefaultScalingType(type);

    FunctionTree<2>::setDefaultOrder(k);
    FunctionTree<2>::setDefaultMaxDepth(depth);
    FunctionTree<2>::setDefaultPrecision(prec);
    FunctionTree<2>::setDefaultScalingType(type);

    FunctionTree<3>::setDefaultOrder(k);
    FunctionTree<3>::setDefaultMaxDepth(depth);
    FunctionTree<3>::setDefaultPrecision(prec);
    FunctionTree<3>::setDefaultScalingType(type);
*/
}

void MREnv::finalizeMRCPP(double t) {
    SET_PRINT_PRECISION(5);
    println(0, endl);
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***                     Exiting MRGrid                   ***");
    println(0, "***                                                      ***");
    println(0, "***               World clock: " << t << "               ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, endl);
}
