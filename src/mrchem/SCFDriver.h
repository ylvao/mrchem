#ifndef SCFDRIVER_H
#define SCFDRIVER_H

#include <vector>
#include <string>
#include <Eigen/Core>

class Getkw;

class Nuclei;
class Molecule;
class OrbitalVector;
class OrbitalOptimizer;
class EnergyOptimizer;
class GroundStateSolver;
class HelmholtzOperatorSet;
class KAIN;

class QMOperator;
class FockOperator;
class CoulombPotential;
class KineticOperator;
class NuclearPotential;
class ExchangePotential;
class XCPotential;
class XCFunctional;

class SCFDriver {
public:
    SCFDriver(Getkw &input);
    virtual ~SCFDriver() { }

    void setup();
    void run();
    void clear();

protected:
    // Top level input
    double rel_prec;

    // World input
    bool center_of_mass;
    std::vector<double> gauge;

    // Run parameters
    bool calc_total_energy;
    bool calc_dipole_moment;

    // Molecule input
    int mol_charge;
    int mol_multiplicity;
    std::vector<std::string> mol_coords;

    // Wavefunction input
    bool wf_restricted;
    std::string wf_method;

    // DFT input
    bool dft_spin;
    double dft_x_fac;
    double dft_cutoff;
    std::vector<double> dft_func_coefs;
    std::vector<std::string> dft_func_names;

    // Ground state input
    std::string scf_start;
    std::string scf_acc;
    int scf_history;
    int scf_max_iter;
    int scf_rotation;
    bool scf_run;
    bool scf_localize;
    bool scf_write_orbitals;
    double scf_orbital_thrs;
    double scf_property_thrs;
    double scf_lambda_thrs;
    std::vector<double> scf_orbital_prec;

    // File input
    std::string file_start_orbitals;
    std::string file_final_orbitals;
    std::string file_basis_set;
    std::string file_dens_mat;
    std::string file_fock_mat;
    std::string file_energy_vec;
    std::string file_mo_mat_a;
    std::string file_mo_mat_b;

    // Gauge origin
    double r_O[3];

    // SCF machinery
    HelmholtzOperatorSet *helmholtz;
    KAIN *scf_kain;

    // Unperturbed quantities
    Molecule *molecule;
    Nuclei *nuclei;
    OrbitalVector *phi;
    KineticOperator *T;
    NuclearPotential *V;
    CoulombPotential *J;
    ExchangePotential *K;
    XCPotential *XC;
    FockOperator *fock;
    Eigen::MatrixXd F;

    OrbitalVector *phi_np1;
    CoulombPotential *J_np1;
    ExchangePotential *K_np1;
    XCPotential *XC_np1;
    FockOperator *fock_np1;
    Eigen::MatrixXd F_np1;

    // XCFun
    XCFunctional *xcfun;

    bool sanityCheck() const;

    bool runGroundState();
    void calcGroundStateProperties();

    void setupInitialGroundState();
    GroundStateSolver *setupInitialGuessSolver();
    OrbitalOptimizer *setupOrbitalOptimizer();
    EnergyOptimizer *setupEnergyOptimizer();

    void setup_np1();
    void clear_np1();

    void printEigenvalues(OrbitalVector &orbs, Eigen::MatrixXd &f_mat);
};

#endif
