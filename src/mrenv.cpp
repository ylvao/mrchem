#include "MRCPP/Printer"

#include <fstream>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

using json = nlohmann::json;
using Printer = mrcpp::Printer;

namespace mrchem {

namespace mrenv {
void init_printer(const json &json_print);
void init_mra(const json &json_mra);
void init_mpi(const json &json_mpi);
void print_header();
void print_footer(double wt);
} // namespace mrenv

json mrenv::fetch_input(int argc, char **argv) {
    const char *infile = nullptr;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    // Read JSON input
    json input;
    std::ifstream ifs(infile, std::ios_base::in);
    ifs >> input;
    ifs.close();

    return input;
}

void mrenv::initialize(const json &input) {
    auto json_print = input.find("printer");
    auto json_mra = input.find("mra");
    auto json_mpi = input.find("mpi");

    if (json_mra == input.end()) MSG_ABORT("Missing MRA input!");

    if (json_print != input.end()) mrenv::init_printer(*json_print);
    if (json_mra != input.end()) mrenv::init_mra(*json_mra);
    if (json_mpi != input.end()) mrenv::init_mpi(*json_mpi);

    mrenv::print_header();
}

void mrenv::init_printer(const json &json_print) {
    // Initialize printing
    auto print_level = json_print["print_level"].get<int>();
    auto print_prec = json_print["print_prec"].get<int>();
    auto mpi_print = json_print["mpi_print"].get<bool>();
    auto fname = json_print["file_name"].get<std::string>();
    if (mpi_print) {
        Printer::init(print_level, mpi::orb_rank, mpi::orb_size, fname.c_str());
    } else {
        Printer::init(print_level, mpi::orb_rank, mpi::orb_size);
    }
    Printer::setPrecision(print_prec);
    Printer::setWidth(60);
}

void mrenv::init_mra(const json &json_mra) {
    // Initialize world box
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
    if (min_scale < mrcpp::MinScale) MSG_ABORT("Root scale too large");
    if (max_scale > mrcpp::MaxScale) MSG_ABORT("Max scale too large");
    if (max_depth > mrcpp::MaxDepth) MSG_ABORT("Max depth too large");

    // Initialize global MRA
    if (btype == "i") {
        mrcpp::InterpolatingBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else if (btype == "l") {
        mrcpp::LegendreBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else {
        MSG_ABORT("Invalid basis type!");
    }
}

void mrenv::init_mpi(const json &json_mpi) {
    auto exact = json_mpi["numerically_exact"].get<bool>();
    auto memory_size = json_mpi["shared_memory_size"].get<int>();
    mpi::numerically_exact = exact;
    mpi::shared_memory_size = memory_size;
}

void mrenv::print_header() {
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

    mrcpp::print::environment(0);
    MRA->print();
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
