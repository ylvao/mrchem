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
#include "parallel.h"

#include "driver.h"
#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "utils/Bank.h"

#include "chemistry/Molecule.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    mpi::initialize();
    const auto json_input = mrenv::fetch_input(argc, argv);

    mrenv::initialize(json_input);
    const auto &json_mol = json_input["molecule"].get<json>();
    const auto &json_guess = json_input["initial_guess"].get<json>();
    const auto &json_scf = json_input["scf_calculation"].get<json>();
    const auto &json_rsps = json_input["rsp_calculations"].get<json>();

    Timer timer;
    Molecule mol;
    if (mpi::is_bankclient) {
        mrcpp::print::memory(2, "memusage before init_molecule");
        driver::init_molecule(json_mol, mol);
        mrcpp::print::memory(2, "memusage before run_guess");
        driver::run_guess(json_guess, mol);
        mrcpp::print::memory(2, "memusage before run_scf");
        if (driver::run_scf(json_scf, mol)) {
            for (const auto &json_rsp : json_rsps) driver::run_rsp(json_rsp, mol);
        }
        driver::print_properties(mol);
        if (mpi::grand_master()) orb_bank.close();
        mpi::barrier(mpi::comm_orb);
        mrenv::finalize(timer.elapsed());
    } else {
        orb_bank.open();
    }

    mpi::finalize();
    return EXIT_SUCCESS;
}
