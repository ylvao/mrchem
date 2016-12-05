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
class LinearResponseSolver;
class HelmholtzOperatorSet;
class KAIN;

class PSOOperator;
class DipoleOperator;
class AngularMomentumOperator;

class QMOperator;
class FockOperator;
class CoulombPotential;
class KineticOperator;
class NuclearPotential;
class ExchangePotential;
class XCPotential;
class XCFunctional;


class ResponseCalculation {
public:
    ResponseCalculation(QMOperator *h, double w, bool im, int d)
        : pert(h), freq(w), imag(im), dir(d) { }

    bool isDynamic() const { if (fabs(this->freq) > MachineZero) return true; return false; }
    bool isImaginary() const { return this->imag; }

    QMOperator *pert;
    double freq;
    bool imag;
    int dir;
};

class ResponseCalculations {
public:
    void push_back(QMOperator *h, double w, bool im, int d) {
        ResponseCalculation rsp_calc(h, w, im, d);
        bool unique = false;
        for (int i = 0; i < this->calculations.size(); i++) {
            const ResponseCalculation &i_calc = calculations[i];
            if (i_calc.pert != rsp_calc.pert) unique = true;
            if (fabs(i_calc.freq - rsp_calc.freq) > MachineZero) unique = true;
            if (i_calc.dir != rsp_calc.dir) unique = true;
        }
        if (unique) {
            this->calculations.push_back(rsp_calc);
        }
    }

    void clear() { this->calculations.clear(); }

    int size() const { return this->calculations.size(); }
    const ResponseCalculation &operator[](int i) const { return this->calculations[i]; }

protected:
    std::vector<ResponseCalculation> calculations;
};

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
    bool calc_scf_energy;
    bool calc_dipole_moment;
    bool calc_quadrupole_moment;
    bool calc_polarizability;
    bool calc_hyperpolarizability;
    bool calc_optical_rotation;
    bool calc_magnetizability;
    bool calc_nmr_shielding;
    bool calc_spin_spin_coupling;
    bool calc_hyperfine_coupling;
    ResponseCalculations rsp_calculations;

    bool pol_velocity;
    bool optrot_velocity;
    string nmr_perturbation;
    string optrot_perturbation;
    std::vector<double> pol_frequency;
    std::vector<double> optrot_frequency;
    std::vector<int> nmr_nucleus_k;
    std::vector<int> sscc_nucleus_k;
    std::vector<int> sscc_nucleus_l;
    std::vector<int> hfcc_nucleus_k;

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

    // Linear response input
    std::string rsp_start;
    int rsp_history;
    int rsp_max_iter;
    bool rsp_run;
    bool rsp_localize;
    bool rsp_write_orbitals;
    double rsp_orbital_thrs;
    double rsp_property_thrs;
    std::vector<int> rsp_directions;
    std::vector<double> rsp_orbital_prec;

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

    // Perturbed quantities
    OrbitalVector *phi_x;
    OrbitalVector *phi_y;
    CoulombPotential *dJ;
    ExchangePotential *dK;
    XCPotential *dXC;
    FockOperator *d_fock;

    // XCFun
    XCFunctional *xcfun;

    // Perturbation operators
    DipoleOperator **h_E;           // dH/dE   [x,y,z]
    AngularMomentumOperator **h_B;  // dH/dB   [x,y,z]
    PSOOperator ***h_M;             // dH/dM_K [x,y,z]

    bool sanityCheck() const;

    bool runGroundState();
    void runLinearResponse(const ResponseCalculation &rsp_calc);

    void calcGroundStateProperties();
    void calcLinearResponseProperties();

    void setupInitialGroundState();
    void setupPerturbedOperators(const ResponseCalculation &rsp_calc);
    void setupPerturbedOrbitals(bool dynamic);

    void clearPerturbedOperators();
    void clearPerturbedOrbitals(bool dynamic);

    GroundStateSolver *setupInitialGuessSolver();
    OrbitalOptimizer *setupOrbitalOptimizer();
    EnergyOptimizer *setupEnergyOptimizer();
    LinearResponseSolver* setupLinearResponseSolver(bool dynamic);

    void setup_np1();
    void clear_np1();

    void printEigenvalues(OrbitalVector &orbs, Eigen::MatrixXd &f_mat);
};

#endif
