#include "MRCPP/Printer"

#include <fstream>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

using json = nlohmann::json;
using Printer = mrcpp::Printer;

namespace mrchem {

void mrenv::initialize(int argc, char **argv) {
    const char *infile = nullptr;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    // Read JSON input
    std::ifstream ifs(infile, std::ios_base::in);
    ifs >> json_input;
    ifs.close();

    // Initialize printing
    auto json_print = json_input["printer"].get<json>();
    auto printlevel = json_print["printlevel"].get<int>();
    auto printprec = json_print["printprec"].get<int>();
    auto teletype = json_print["teletype"].get<bool>();
    auto fname = json_print["filename"].get<std::string>();
    if (teletype) {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size, fname.c_str());
    } else {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size);
    }
    Printer::setPrecision(printprec);

    // Initialize global MPI parameters
    auto json_mpi = json_input["mpi"].get<json>();
    mpi::numerically_exact = json_mpi["numerically_exact"].get<bool>();
    mpi::share_nuc_pot = json_mpi["share_nuclear_potential"].get<bool>();
    mpi::share_coul_dens = json_mpi["share_coulomb_density"].get<bool>();
    mpi::share_coul_pot = json_mpi["share_coulomb_potential"].get<bool>();
    mpi::share_xc_dens = json_mpi["share_xc_density"].get<bool>();
    mpi::share_xc_pot = json_mpi["share_xc_potential"].get<bool>();
    mpi::shared_memory_size = json_mpi["shared_memory_size"].get<int>();

    // Initialize world box
    auto json_mra = json_input["mra"].get<json>();
    auto min_scale = json_mra["min_scale"].get<int>();
    auto max_scale = json_mra["max_scale"].get<int>();
    auto corner = json_mra["corner"].get<std::array<int, 3>>();
    auto boxes = json_mra["boxes"].get<std::array<int, 3>>();
    auto sfac = json_mra["scaling_factor"].get<std::array<double, 3>>();
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes, sfac);

    // Initialize scaling basis
    auto order = json_mra["order"].get<int>();
    auto btype = json_mra["basis_type"].get<std::string>();

    auto max_depth = max_scale - min_scale;
    if (min_scale < mrcpp::MinScale) MSG_FATAL("Root scale too large");
    if (max_scale > mrcpp::MaxScale) MSG_FATAL("Max scale too large");
    if (max_depth > mrcpp::MaxDepth) MSG_FATAL("Max depth too large");

    // Initialize global MRA
    if (btype == "i") {
        mrcpp::InterpolatingBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else if (btype == "l") {
        mrcpp::LegendreBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else {
        MSG_FATAL("Invalid basis type!");
    }

    println(0, std::endl << std::endl);
    println(0, "************************************************************");
    println(0, "***     __  __ ____   ____ _                             ***");
    println(0, "***    |  \\/  |  _ \\ / ___| |__   ___ _ __ ___           ***");
    println(0, "***    | |\\/| | |_) | |   | '_ \\ / _ \\ '_ ` _ \\          ***");
    println(0, "***    | |  | |  _ <| |___| | | |  __/ | | | | |         ***");
    println(0, "***    |_|  |_|_| \\_\\\\____|_| |_|\\___|_| |_| |_|         ***");
    println(0, "***                                                      ***");
    println(0, "***    VERSION " << PROGRAM_VERSION << " (rev. " << GIT_REVISION << ")                     ***");
    println(0, "***                                                      ***");
    println(0, "***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0, "***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0, "***    Peter Wind       <peter.wind@uit.no>              ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, std::endl);

    if (mpi::orb_size > 1 or omp::n_threads > 1) {
        println(0, "+++ Parallel execution: ");
        println(0, "  MPI hosts available     : " << mpi::orb_size);
        println(0, "  Threads/host            : " << omp::n_threads);
        println(0, "  Total used CPUs         : " << mpi::orb_size * omp::n_threads);
        println(0, "");
    } else {
        println(0, "+++ Serial execution" << std::endl);
    }

    Printer::printEnvironment();

    Printer::printHeader(0, "JSON input");
    println(0, std::setw(2) << json_input);
    Printer::printSeparator(0, '=', 2);

    auto oldprec = Printer::setPrecision(6);
    MRA->print();
    Printer::setPrecision(oldprec);
}

void mrenv::finalize(double wt) {
    // Delete global MRA
    if (MRA != nullptr) delete MRA;
    MRA = nullptr;

    auto oldprec = Printer::setPrecision(6);
    println(0, std::endl);
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***                    Exiting MRChem                    ***");
    println(0, "***                                                      ***");
    println(0, "***               Wall time: " << wt << "                ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, std::endl);
    Printer::setPrecision(oldprec);
}

} // namespace mrchem
