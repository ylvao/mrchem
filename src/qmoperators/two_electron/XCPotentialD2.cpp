#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "XCPotential.h"
#include "XCPotentialD2.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

int XCPotentialD2::getPotentialIndex(int orbitalSpin, int densitySpin);

    /** @brief Constructor
 *
 * @param[in] F XCFunctional pointer
 * @param[in] Phi Vector of orbitals
 *
 * Based on the order and spin the correct nr. of potential functions is determined.
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 */
XCPotentialD2::XCPotentialD2(mrdft::XCFunctional *F, OrbitalVector *Phi)
        : XCPotential(F)
          orbitals(Phi),
          orbitals_x(X),
          orbitals_y(Y) {
}

/** @brief Prepare the operator for application
 * 
 * @param[in] prec Apply precision
 *
 * Sequence of steps required to compute the XC potentials:
 *
 * 1) Compute density
 * 2) Setup xcfun input functions (gradients etc.)
 * 3) Evaluate xcfun
 * 4) Compute XC potential(s) from xcfun output
 *
 */
void XCPotentialD2::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity(prec);
    setupPotential(prec);
}

/** @brief Clears all data in the XCPotentialD2 object */
void XCPotentialD2::clear() {
    this->energy = 0.0;
    mrcpp::clear(this->potentials, true);
    clearApplyPrec();
}

/** @brief Compute electron density
 *
 * The density is computed on the grid provided by the MRDFT module. The grid
 * is kept as is, e.i. no additional refinement at this point, since the grid
 * size is determined inside the module.
 */
void XCPotentialD2::setupDensity() {
    setupGroundStateDensity();
    setupPerturbedDensity();
}

void XCPotentialD2::setupGroundStateDensity() {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        FunctionTree<3> &tmp_a = getDensity(DENSITY::Alpha);
        Density rho_a;
        rho_a.setReal(&tmp_a);
        density::compute(-1.0, rho_a, Phi, DENSITY::Alpha);
        time_a.stop();
        Printer::printTree(0, "XC alpha density", rho_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        FunctionTree<3> &tmp_b = getDensity(DENSITY::Beta);
        Density rho_b;
        rho_b.setReal(&tmp_b);
        density::compute(-1.0, rho_b, Phi, DENSITY::Beta);
        time_b.stop();
        Printer::printTree(0, "XC beta density", rho_b.getNNodes(), time_b.getWallTime());
    } else {
        Timer time_t;
        FunctionTree<3> &tmp_t = getDensity(DENSITY::Total);
        Density rho_t;
        rho_t.setReal(&tmp_t);
        density::compute(-1.0, rho_t, Phi, DENSITY::Total);
        time_t.stop();
        Printer::printTree(0, "XC total density", rho_t.getNNodes(), time_t.getWallTime());
    }
}

void XCPotentialD2::setupPerturbedDensity() {
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        FunctionTree<3> &tmp_a = getDensity(DENSITY::Alpha);
        pertDensity_a->allocReal();
        mrcpp::copy_grid(pertDensity_a->real(), tmp_a);
        density::compute(-1.0, pertDensity_a, Phi, X, Y, DENSITY::Alpha);
        time_a.stop();
        Printer::printTree(0, "XC alpha density", pertDensity_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        FunctionTree<3> &tmp_b = getDensity(DENSITY::Beta);
        pertDensity_b->allocReal();
        mrcpp::copy_grid(pertDensity_b->real(), tmp_b);
        density::compute(-1.0, pertDensity_b, Phi, X, Y, DENSITY::Beta);
        time_b.stop();
        Printer::printTree(0, "XC beta density", pertDensity_b.getNNodes(), time_b.getWallTime());

    } else {
        Timer time_t;
        FunctionTree<3> &tmp_t = getDensity(DENSITY::Total);
        pertDensity_t->allocReal();
        mrcpp::copy_grid(pertDensity_t->real(), tmp_t);
        density::compute(-1.0, pertDensity_t, Phi, X, Y, DENSITY::Total);
        time_t.stop();
        Printer::printTree(0, "XC total density", pertDensity_t.getNNodes(), time_t.getWallTime());
    }
}

/** @brief Compute XC potential(s)
 *
 * @param prec Precision used in refinement of density grid
 *
 * This will invoke a sequence of steps in the XCFunctional to compute the final
 * XC potential(s) that define this operator. Assuming the density has alredy been
 * computed:
 *
 * 1) Setup xcfun input functions (gradients etc.)
 * 2) Evaluate xcfun
 * 3) Compute XC energy by integrating energy density
 * 4) Compute XC potential(s) from xcfun output functions
 * 5) Remove excess grid nodes based on precision
 * 6) Add extra grid nodes based on precision
 * 7) Clear internal functions in XCFunctional (density grid is kept)
 *
 */
void XCPotentialD2::setupPotential(double prec) {
    if (this->functional == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (not this->functional->hasDensity()) MSG_ERROR("XC density not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    int inpNodes = this->functional->getNNodes();
    int inpPoints = this->functional->getNPoints();

    this->functional->setup();
    this->functional->evaluate();
    //    this->energy = this->functional->calcEnergy();
    this->potentials = this->functional->calcPotential();
    this->functional->pruneGrid(prec);
    this->functional->refineGrid(prec);
    this->functional->clear();

    int newNodes = this->functional->getNNodes() - inpNodes;
    int newPoints = this->functional->getNPoints() - inpPoints;

    println(0, " XC grid size   " << std::setw(26) << inpNodes << std::setw(17) << inpPoints);
    println(0, " XC grid change " << std::setw(26) << newNodes << std::setw(17) << newPoints);
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
 */
FunctionTree<3>& XCPotentialD2::getPotential(int orbitalSpin, int densitySpin) {
    int pot_idx = XCPotentialD2::getPotentialIndex(int orbitalSpin, int densitySpin);
    return mrcpp::get_func(this->potentials, pot_idx);
}

int XCPotentialD2::getPotentialIndex(int orbitalSpin, int densitySpin) {

    int spinFunctional = this->functional->isSpinSeparated()?:1:0;
    
    int functional_case;
    functional_case += spinFunctional;       // 0  1
    functional_case += orbitalSpin   << 2;   // 0  2  4  6
    functional_case += densitySpin   << 4;   // 0  8 16 24

    //OS Paired Alpha Beta
    //DS Total  Spin  Alpha Beta
    
    // SF    OS             DS          case   Index 
    //  0     0 (paired)     0 (total)     0    0    
    //  1     0 (paired)     0 (total)     1    not implemented
    //  0     1 (alpha )     0 (total)     2    not implemented 
    //  1     1 (alpha )     0 (total)     3    not implemented 
    //  0     2 (beta  )     0 (total)     4    not implemented 
    //  1     2 (beta  )     0 (total)     5    not implemented 
    //  0     3 (unused)     0 (total)     6    not implemented 
    //  1     3 (unused)     0 (total)     7    not implemented 
    //  0     0 (paired)     1 (spin )     8    not implemented 
    //  1     0 (paired)     1 (spin )     9    not implemented 
    //  0     1 (alpha )     1 (spin )    10    not implemented 
    //  1     1 (alpha )     1 (spin )    11    not implemented 
    //  0     2 (beta  )     1 (spin )    12    not implemented 
    //  1     2 (beta  )     1 (spin )    13    not implemented 
    //  0     3 (unused)     1 (spin )    14    not implemented  
    //  1     3 (unused)     1 (spin )    15    not implemented  
    //  0     0 (paired)     2 (alpha)    16    not implemented  
    //  1     0 (paired)     2 (alpha)    17    not implemented  
    //  0     1 (alpha )     2 (alpha)    18    not implemented  
    //  1     1 (alpha )     2 (alpha)    19    not implemented  
    //  0     2 (beta  )     2 (alpha)    20    not implemented  
    //  1     2 (beta  )     2 (alpha)    21    not implemented  
    //  0     3 (unused)     2 (alpha)    22    not implemented  
    //  1     3 (unused)     2 (alpha)    23    not implemented  
    //  0     0 (paired)     3 (beta )    24    not implemented  
    //  1     0 (paired)     3 (beta )    25    not implemented  
    //  0     1 (alpha )     3 (beta )    26    not implemented  
    //  1     1 (alpha )     3 (beta )    27    not implemented  
    //  0     2 (beta  )     3 (beta )    28    not implemented  
    //  1     2 (beta  )     3 (beta )    29    not implemented  
    //  0     3 (unused)     3 (beta )    30    not implemented  
    //  1     3 (unused)     3 (beta )    31    not implemented  
    switch(functional_case) {
    case(0) potential_index = 0;    break;
    default: MSG_FATAL("Not implemented: ABORT");
    }
}

/** @brief XCPotentialD2 application
 *
 * @param[in] phi Orbital to which the potential is applied
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 */
Orbital XCPotentialD2::apply(Orbital phi) {
    if (this->hasImag()) MSG_ERROR("Imaginary part of XC potential non-zero");

    FunctionTree<3> &V = getPotential(phi.spin());
    this->setReal(&V);
    Orbital Vphi = QMPotential::apply(phi); 
    this->setReal(0);

    return Vphi;
}

} //namespace mrchem
