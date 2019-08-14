#include "MRCPP/Printer"

#include <fstream>

#include "XCFun/xcfun.h"

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "utils/print_utils.h"

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
    auto print_width = json_print["print_width"].get<int>();
    auto print_mpi = json_print["print_mpi"].get<bool>();
    auto fname = json_print["file_name"].get<std::string>();
    if (print_mpi) {
        Printer::init(print_level, mpi::world_rank, mpi::world_size, fname.c_str());
    } else {
        Printer::init(print_level, mpi::world_rank, mpi::world_size);
    }
    Printer::setPrecision(print_prec);
    Printer::setWidth(print_width);
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
    auto bank_size_json = json_mpi["bank_size"].get<int>();
    mpi::numerically_exact = exact;
    mpi::shared_memory_size = memory_size;
    mpi::bank_size = bank_size_json;
}

void mrenv::print_header() {
    auto pwidth = Printer::getWidth();
    auto txt_width = 45;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');
    std::stringstream o_ver, o_rev;
    o_ver << "VERSION " << std::setw(8) << PROGRAM_VERSION;
    o_rev << "(rev. " << std::setw(10) << GIT_REVISION << ")";

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << " __  __ ____   ____ _                        " << post_str);
    println(0, pre_str << "|  \\/  |  _ \\ / ___| |__   ___ _ __ ___      " << post_str);
    println(0, pre_str << "| |\\/| | |_) | |   | '_ \\ / _ \\ '_ ` _ \\     " << post_str);
    println(0, pre_str << "| |  | |  _ <| |___| | | |  __/ | | | | |    " << post_str);
    println(0, pre_str << "|_|  |_|_| \\_\\\\____|_| |_|\\___|_| |_| |_|    " << post_str);
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << o_ver.str() << "        " << o_rev.str() << "    " << post_str);
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "Stig Rune Jensen   <stig.r.jensen@uit.no>    " << post_str);
    println(0, pre_str << "Luca Frediani      <luca.frediani@uit.no>    " << post_str);
    println(0, pre_str << "Peter Wind         <peter.wind@uit.no>       " << post_str);
    println(0, pre_str << "                                             " << post_str);
    mrcpp::print::separator(0, '*', 1);
    mrcpp::print::separator(0, '-', 1);
    print_utils::scalar(0, "MPI processes", mpi::world_size, "", 0, false);
    print_utils::scalar(0, "of which used as bank", mpi::bank_size, "", 0, false);
    print_utils::scalar(0, "of which used for orb", mpi::orb_size, "", 0, false);
    print_utils::scalar(0, "OpenMP threads", omp::n_threads, "", 0, false);
    print_utils::scalar(0, "Total cores", mpi::orb_size * omp::n_threads, "", 0, false);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '-', 1);
    printout(0, xcfun_splash());
    mrcpp::print::environment(0);
    MRA->print();
}

void mrenv::finalize(double wt) {
    // Delete global MRA
    if (MRA != nullptr) delete MRA;
    MRA = nullptr;

    auto pwidth = Printer::getWidth();
    auto txt_width = 45;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    auto hr = static_cast<int>(wt / 3600.0);
    auto min = static_cast<int>(std::fmod(wt, 3600.0) / 60.0);
    auto sec = static_cast<int>(std::fmod(wt, 60.0));
    std::stringstream o_time;
    o_time << "Wall time : " << std::setw(2) << hr << "h" << std::setw(3) << min << "m" << std::setw(3) << sec << "s";

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "                Exiting MRChem               " << post_str);
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "           " << o_time.str() << "           " << post_str);
    println(0, pre_str << "                                             " << post_str);
    mrcpp::print::separator(0, '*');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
}

} // namespace mrchem
