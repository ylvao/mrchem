#include "surface_forces/SurfaceForce.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/NablaOperator.h"
#include <vector>
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PhysicalConstants.h"
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

#include "mrdft/Factory-orig.h"
#include "mrdft/MRDFT.h"
#include "mrdft/xc_utils.h"
#include "qmoperators/two_electron/XCOperator.h"
#include "mrdft/Functional.h"

#include "qmoperators/one_electron/HessianOperator.h"
#include "tensor/RankOneOperator.h"

#include "surface_forces/lebedev.h"
#include "surface_forces/xcStress.h"
#include <string>
#include <iostream>
#include <filesystem>

#include <MRCPP/Timer>

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace Eigen;
using namespace mrchem;
using nlohmann::json;

namespace surface_force {

/**
 * Calculates the electric field due to the nuclei.
 * @param nucPos The positions of the nuclei. Shape (nNuc, 3).
 * @param nucCharge The charges of the nuclei. Shape (nNuc,).
 * @param nucSmoothing The smoothing parameter for the nuclei. Shape (nNuc,).
 * @param gridPos The positions of the grid points where the field should be evaluated. Shape (nGrid, 3).
*/
MatrixXd nuclearEfield(const MatrixXd &nucPos, const VectorXd &nucCharge, const VectorXd &nucSmoothing, const MatrixXd gridPos) {
    int nGrid = gridPos.rows();
    int nNuc = nucPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    Vector3d r_vect;
    double r;
    double r2, r3;
    double temp;
    double c, q;
    double c2, c3;
    double c3_times_sqrt_pi_times_three;
    double sqrt_pi = std::sqrt(M_PI);
    for (int i = 0; i < nNuc; ++i) {
        c = nucSmoothing(i);
        q = nucCharge(i);
        c2 = c * c;
        c3 = c2 * c;
        c3_times_sqrt_pi_times_three = 3. * std::sqrt(M_PI) * c3;
        for (int j = 0; j < nGrid; j++)
        {
            r_vect = nucPos.row(i) - gridPos.row(j);
            r = r_vect.norm();
            r2 = r * r;
            r3 = r2 * r;
            Efield.row(j) += q * r_vect / r3;
            // This would compute the electric field for the finite point like nucleus. It does not improve accuracy since the integration speres
            // are large enough to not be affected by the finite point like nucleus. It is kept here for reference.
            // Efield.row(j) += q * r_vect / r * (std::erf(r / c) / r2 - 2 * std::exp(-r2/c2)/(sqrt_pi * c * r) + 2.0 * r * std::exp(-r2 / c2) / c3_times_sqrt_pi_times_three 
            //     + 128.0 * r * std::exp(-4.0 * r2 / c2) / c3_times_sqrt_pi_times_three);
        }
    }
    return Efield;
}

/**
 * @brief Calculates the coulomb potential due to the given density.
*/
mrcpp::ComplexFunction calcPotential(Density &rho, mrcpp::PoissonOperator &poisson, double prec) {
    mrcpp::ComplexFunction V(false);
    V.alloc(mrchem::NUMBER::Real);
    mrcpp::apply(prec, V.real(), poisson, rho.real());
    return V;
}

/**
 * @brief Calculates the electric field due to the electrons.
 * @param negEfield The negative gradient of the electric field. Shape (3,).
 * @param gridPos The positions of the grid points where the field should be evaluated. Shape (nGrid, 3).
*/
MatrixXd electronicEfield(mrchem::OrbitalVector &negEfield, const MatrixXd &gridPos) {
    int nGrid = gridPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    for (int i = 0; i < nGrid; i++)
    {
        std::array<double, 3> pos = {gridPos(i, 0), gridPos(i, 1), gridPos(i, 2)};
        Efield(i, 0) = -negEfield[0].real().evalf(pos);
        Efield(i, 1) = -negEfield[1].real().evalf(pos);
        Efield(i, 2) = -negEfield[2].real().evalf(pos);
    }
    return Efield;
}

/**
 * Calculates the Maxwell stress tensor for the given molecule.
 * @param mol The molecule for which to calculate the stress tensor.
 * @param negEfield Negative electric field (gradient of potential)
 * @param gridPos The positions of the grid points where the field should be evaluated. Shape (nGrid, 3).
 * @param prec The precision value used in the calculation.
*/
std::vector<Eigen::Matrix3d> maxwellStress(const Molecule &mol, mrchem::OrbitalVector &negEfield, const MatrixXd &gridPos, double prec){
    int nGrid = gridPos.rows();
    int nNuc = mol.getNNuclei();

    Eigen::MatrixXd nucPos(nNuc, 3);
    Eigen::VectorXd nucCharge(nNuc);
    Eigen::VectorXd nucSmoothing(nNuc);
    for (int i = 0; i < nNuc; i++)
    {
        std::array<double, 3> coord = mol.getNuclei()[i].getCoord();
        nucPos(i, 0) = coord[0];
        nucPos(i, 1) = coord[1];
        nucPos(i, 2) = coord[2];
        nucCharge(i) = mol.getNuclei()[i].getCharge();
        double tmp = 0.00435 * prec / std::pow(nucCharge(i), 5.0);
        nucSmoothing(i) =  std::cbrt(tmp);
    }

    MatrixXd Efield = electronicEfield(negEfield, gridPos) + nuclearEfield(nucPos, nucCharge, nucSmoothing, gridPos);

    std::vector<Eigen::Matrix3d> stress(nGrid);
    for (int i = 0; i < nGrid; i++) {
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress[i](i1, i2) = Efield(i, i1) * Efield(i, i2);
            }
        }
        for (int i1 = 0; i1 < 3; i1++){
            stress[i](i1, i1) = stress[i](i1, i1) - 0.5 * (Efield(i, 0) * Efield(i, 0) + Efield(i, 1) * Efield(i, 1) + Efield(i, 2) * Efield(i, 2));
        }
        for (int i1 = 0; i1 < 3; i1++){
            for (int i2 = 0; i2 < 3; i2++){
                stress[i](i1, i2) *= 1.0 / (4 * M_PI);
            }
        }
    }
    return stress;
}

/**
 * @brief Calculates the kinetic stress tensor for the given molecule. See the function description for the formula.
*/
std::vector<Matrix3d> kineticStress(const Molecule &mol, OrbitalVector &Phi, std::vector<std::vector<mrchem::Orbital>> &nablaPhi
        , std::vector<Orbital> &hessRho, double prec, const MatrixXd &gridPos){

    // original formula for kinetic stress:
    // sigma_ij = 0.5 \sum_k phi_k del_i del_j phi_k - (del_i phi_k) (del_j phi_k)
    // using the product rule for the first term:
    // sigma_ij = 0.5 * del_i del_j rho - 2 * \sum_k (del_i phi_k) (del_j phi_k)
    // That way, only second derivatives of density is needed which speeds things up. quite a lot.

    int nGrid = gridPos.rows();
    int nOrbs = Phi.size();

    double StressArray[nGrid][3][3];

    double orbVal;
    std::vector<Matrix3d> stress(nGrid);

    Eigen::MatrixXd voigtStress = Eigen::MatrixXd::Zero(nGrid, 6);

    std::array<double, 3> pos;
    double n1, n2, n3;
    int occ;
    for (int iOrb = 0; iOrb < Phi.size(); iOrb++) {
        occ = Phi[iOrb].occ();
    
        for (int i = 0; i < nGrid; i++) {
            if (mrcpp::mpi::my_orb(iOrb)) {  
            pos[0] = gridPos(i, 0);
            pos[1] = gridPos(i, 1);
            pos[2] = gridPos(i, 2);
            n1 = nablaPhi[iOrb][0].real().evalf(pos);
            n2 = nablaPhi[iOrb][1].real().evalf(pos);
            n3 = nablaPhi[iOrb][2].real().evalf(pos);
            voigtStress(i, 0) -= occ * n1 * n1;
            voigtStress(i, 1) -= occ * n2 * n2;
            voigtStress(i, 2) -= occ * n3 * n3;
            voigtStress(i, 5) -= occ * n1 * n2;
            voigtStress(i, 4) -= occ * n1 * n3;
            voigtStress(i, 3) -= occ * n2 * n3;
            }

        }
    }
    mrcpp::mpi::allreduce_matrix(voigtStress, mrcpp::mpi::comm_wrk);
    for (int i = 0; i < nGrid; i++) {
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        voigtStress(i, 0) += 0.25 * hessRho[0].real().evalf(pos);
        voigtStress(i, 1) += 0.25 * hessRho[1].real().evalf(pos);
        voigtStress(i, 2) += 0.25 * hessRho[2].real().evalf(pos);
        voigtStress(i, 3) += 0.25 * hessRho[3].real().evalf(pos);
        voigtStress(i, 4) += 0.25 * hessRho[4].real().evalf(pos);
        voigtStress(i, 5) += 0.25 * hessRho[5].real().evalf(pos);
    }

    for (int i = 0; i < nGrid; i++) {
        stress[i] << voigtStress(i, 0), voigtStress(i, 5), voigtStress(i, 4),
                     voigtStress(i, 5), voigtStress(i, 1), voigtStress(i, 3),
                     voigtStress(i, 4), voigtStress(i, 3), voigtStress(i, 2);
    }
    return stress;
}

/**
 * Calculates the distance to the nearest neighbor for each point in the given position matrix.
 *
 * @param pos The matrix containing the positions of the points. Shape (nPoints, 3).
 * @return A vector containing the distances to the nearest neighbor for each point.
 */
VectorXd distanceToNearestNeighbour(MatrixXd pos){
    int n = pos.rows();
    VectorXd dist(n);
    double temp;
    if (n == 1){
        dist(0) = 1.0;
    } else {
        for (int i = 0; i < n; i++){
            dist(i) = (pos.row(i) - pos.row((i + 1) % n)).norm();
            for (int j = 0; j < n; j++){
                if (i != j){
                    temp = (pos.row(i) - pos.row(j)).norm();
                    if (temp < dist(i)){
                        dist(i) = temp;
                    }
                }
            }
        }
    }
    return dist;
}

/**
 * Calculates the forces using surface integrals for a given molecule and orbital vector.
 *
 * @param mol The molecule for which to calculate the forces.
 * @param Phi The orbital vector obtained from the SCF calculation.
 * @param prec The precision value used in the calculation.
 * @param json_fock The JSON object containing the Fock matrix settings.
 * @return The matrix of forces, shape (nAtoms, 3).
 */
Eigen::MatrixXd surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const json &json_fock
        , std::string leb_prec, double radius_factor) {

    if (radius_factor > 0.95 && radius_factor < 0.05){
        MSG_ABORT("Invalid value of radius_factor")
    }

    // setup density
    mrchem::Density rho(false);
    mrchem::density::compute(prec, rho, Phi, DensityType::Total);

    std::shared_ptr<OrbitalVector> phi_p = std::make_shared<OrbitalVector>(Phi);

    // setup operators and potentials
    auto poisson_op = mrcpp::PoissonOperator(*mrchem::MRA, prec);
    int derivOrder = 1;
    auto mrcd = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder);
    mrchem::NablaOperator nabla(mrcd);
    nabla.setup(prec);
    double abs_prec = prec / mrchem::orbital::get_electron_number(Phi);
    mrcpp::ComplexFunction pot = calcPotential(rho, poisson_op, abs_prec);
    mrchem::OrbitalVector negEfield = nabla(pot);

    // set up operators for kinetic stress:
    int derivOrder1 = 1;
    auto D1 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder1);
    int derivOrder2 = 2;
    auto D2 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder2);
    mrchem::HessianOperator hess(D1, D2, prec);
    hess.setup(prec);

    std::vector<std::vector<Orbital>> nablaPhi(Phi.size());
    std::vector<Orbital> hessRho = hess(rho);
    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_orb(i)) {
            nablaPhi[i] = nabla(Phi[i]);
        }
    }
    // setup xc stuff:
    int order = 0;
    bool shared_memory = json_fock["xc_operator"]["shared_memory"];
    auto json_xcfunc = json_fock["xc_operator"]["xc_functional"];
    bool xc_spin = json_xcfunc["spin"];
    auto xc_cutoff = json_xcfunc["cutoff"];
    auto xc_funcs = json_xcfunc["functionals"];
    auto xc_order = order + 1;
    auto funcVectorShared = std::make_shared<mrcpp::MPI_FuncVector>(Phi);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(xc_spin);
    xc_factory.setOrder(xc_order);
    xc_factory.setDensityCutoff(xc_cutoff);
    for (const auto &f : xc_funcs) {
        auto name = f["name"];
        auto coef = f["coef"];
        xc_factory.setFunctional(name, coef);
    }
    std::unique_ptr<mrdft::MRDFT> mrdft_p = xc_factory.build();

    mrdft::Factory xc_factory2(*MRA);
    xc_factory2.setSpin(xc_spin);
    xc_factory2.setOrder(xc_order);
    xc_factory2.setDensityCutoff(xc_cutoff);
    for (const auto &f : xc_funcs) {
        auto name = f["name"];
        auto coef = f["coef"];
        xc_factory2.setFunctional(name, coef);
    }
    std::unique_ptr<mrdft::MRDFT> mrdft_p_for_operator = xc_factory2.build();

    std::shared_ptr<XCOperator> xc_op = std::make_shared<XCOperator>(mrdft_p_for_operator, phi_p, false);
    xc_op->setup(prec);

    std::shared_ptr<mrcpp::FunctionTreeVector<3>> xc_pot_vector = xc_op->getPotential()->getPotentialVector();

    int numAtoms = mol.getNNuclei();

    std::array<double, 3> coord;

    Eigen::MatrixXd posmatrix = Eigen::MatrixXd::Zero(numAtoms, 3);
    for (int i = 0; i < numAtoms; i++) {
        coord = mol.getNuclei()[i].getCoord();
        posmatrix(i, 0) = coord[0];
        posmatrix(i, 1) = coord[1];
        posmatrix(i, 2) = coord[2];
    }
    VectorXd dist = distanceToNearestNeighbour(posmatrix);

    int nLebPoints = 0;
    int nTinyPoints = 1;
    if (leb_prec == "low") {
        nLebPoints = 194;
    }
    else if (leb_prec == "medium"){
        nLebPoints = 434;
    }
    else if (leb_prec == "high") {
        nLebPoints = 770;
    }
    else {
        MSG_ABORT("Invalid lebedev precision");
    }

    double radius;
    Vector3d center;
    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(numAtoms, 3);

    for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
        radius = dist(iAtom) * radius_factor;
        coord = mol.getNuclei()[iAtom].getCoord();
        center << coord[0], coord[1], coord[2];

        LebedevIntegrator integrator(nLebPoints, radius, center);
        MatrixXd gridPos = integrator.getPoints();
        VectorXd weights = integrator.getWeights();
        MatrixXd normals = integrator.getNormals();
        std::vector<Matrix3d> xcStress = getXCStress(mrdft_p, *xc_pot_vector, std::make_shared<mrchem::OrbitalVector>(Phi), 
            std::make_shared<mrchem::NablaOperator>(nabla), gridPos, xc_spin, prec);
        
        
        std::vector<Matrix3d> kstress = kineticStress(mol, Phi, nablaPhi, hessRho, prec, gridPos);
        std::vector<Matrix3d> mstress = maxwellStress(mol, negEfield, gridPos, prec);
        std::vector<Matrix3d> stress(integrator.n);
        for (int i = 0; i < integrator.n; i++){
            stress[i] = xcStress[i] + kstress[i] + mstress[i];
            forces.row(iAtom) -= stress[i] * normals.row(i).transpose() * weights(i);
        }
    }

    hess.clear();
    nabla.clear();
    xc_op->clear();
    return forces;
}

} // namespace surface_force
