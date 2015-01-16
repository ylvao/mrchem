#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/timer.hpp>

#include "mrgrid.h"
#include "parallel.h"
#include "config.h"

#include "Getkw.h"
#include "TelePrompter.h"

#include "MREnv.h"
#include "BoundingBox.h"
#include "QuadratureCache.h"

void MREnv::initializeMRCPP(int argc, char **argv, const char *fname) {
    omp_set_dynamic(0);
    mpi::communicator world;

    extern mpi::communicator orbital_group;

    int nThreads = omp_get_max_threads();
    Eigen::setNbThreads(1); 

    const char *infile = 0;
    if (argc == 1) {
	infile = "STDIN";
    } else if (argc == 2) {
	infile = argv[1];
    } else {
	MSG_ERROR("Ivalid number of arguments!");
    }
    //Input = Getkw(infile, false, true);

    Debug = false;
    int printLevel = 0;//Input.get<int>("printlevel");
    bool teletype = 0;//Input.get<bool>("teletype");

    if (fname != 0) {
	TelePrompter::init(printLevel, teletype, fname);
    } else {
	TelePrompter::init(printLevel, teletype, "MRGRID");
    }

    if (printLevel != 0) {
	Debug = true;
    }

    println(0, endl << endl);
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***    MRGrid " << PROJECT_VERSION << " (rev. " <<
	    GIT_REVISION <<                     ")                       ***");
    println(0, "***                                                      ***");
    println(0, "***    Jonas Juselius   <jonas.juselius@uit.no>          ***");
    println(0, "***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
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
    //int order = Input.get<int>("order");
    //int max_depth = Input.get<int>("max_depth");
    //double rel_prec = Input.get<double>("rel_prec");
    //string wlet = Input.get<string>("wavelet");

    //int rootScale = Input.get<int>("World.scale");
    //const vector<int> &nbox = Input.getIntVec("World.boxes");
    //const vector<int> &transl = Input.getIntVec("World.translation");
    //const vector<double> &origin = Input.getDblVec("World.origin");

    //BoundingBox<1>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    //BoundingBox<2>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    //BoundingBox<3>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());

    //const BoundingBox<3> &worldbox = BoundingBox<3>::getWorldBox();
    //println(0, worldbox);

    //int polytype;
    //if (wlet == "I") {
	//polytype = Interpol;
    //} else {
	//polytype = Legendre;
    //}
    //println(0, "*Default parameters:");
    //println(0, "  Debug level  :     " << printLevel);
    //println(0, "  Default order:     " << order);
    //println(0, "  Default max depth: " << max_depth);
    //println(0, "  Default precision: " << rel_prec);
    //printout(0, "  Default polynomial type: ");
    //if (polytype == Interpol) println(1, "Interpolating");
    //if (polytype == Legendre) println(1, "Legendre");
    //println(0, endl);

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
