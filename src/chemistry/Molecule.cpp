/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <fstream>

#include "Molecule.h"
#include "Nucleus.h"

#include "qmfunctions/orbital_utils.h"

using mrcpp::Coord;
using mrcpp::Printer;

namespace mrchem {

/** @brief Constructor
 *
 * @param c: total charge
 * @param m: spin multiplicity
 */
Molecule::Molecule(int c, int m)
        : charge(c)
        , multiplicity(m) {}

/** @brief Constructor
 *
 * @param coord_file: xyz file with nuclear coordinates
 * @param c: total charge
 * @param m: spin multiplicity
 *
 * Nuclei are copied, all properties are uninitialized at this point.
 */
Molecule::Molecule(const std::string &coord_file, int c, int m)
        : charge(c)
        , multiplicity(m) {
    readCoordinateFile(coord_file);
}

/** @brief Constructor
 *
 * @param coord_file: list of stings with nuclear coordinates
 * @param c: total charge
 * @param m: spin multiplicity
 *
 * Nuclei are copied, all properties are uninitialized at this point.
 */
Molecule::Molecule(const std::vector<std::string> &coord_str, int c, int m)
        : charge(c)
        , multiplicity(m) {
    readCoordinateString(coord_str);
}

void Molecule::initNuclearProperties(int nNucs) {
    for (auto k = 0; k < nNucs; k++) {
        this->nmr.push_back(nullptr);
        this->hfcc.push_back(nullptr);
        this->sscc.push_back(std::vector<std::unique_ptr<SpinSpinCoupling>>{});
        for (auto &sscc_k : this->sscc) sscc_k.push_back(nullptr);
    }
}

/** @brief Initialize property SCFEnergy */
void Molecule::initSCFEnergy() {
    if (this->energy != nullptr) MSG_WARN("SCFEnergy already initialized");
    this->energy = std::make_unique<SCFEnergy>();
}

/** @brief Initialize property DipoleMoment */
void Molecule::initDipoleMoment() {
    if (this->dipole != nullptr) MSG_WARN("Dipole moment already initialized");
    this->dipole = std::make_unique<DipoleMoment>();
}

/** @brief Initialize property QuadrupoleMoment */
void Molecule::initQuadrupoleMoment() {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (this->quadrupole != nullptr) MSG_WARN("Quadrupole moment already initialized");
    this->quadrupole = new QuadrupoleMoment();
    */
}

/** @brief Initialize property GeometryDerivatives */
void Molecule::initGeometryDerivatives() {
    if (this->geomderiv != nullptr) MSG_WARN("Geometry derivatives already initialized");
    this->geomderiv = std::make_unique<GeometryDerivatives>(this->nuclei.size());
}

/** @brief Initialize property Magnetizability */
void Molecule::initMagnetizability() {
    if (this->magnetizability != nullptr) MSG_WARN("Magnetizability already initialized");
    this->magnetizability = std::make_unique<Magnetizability>();
}

/** @brief Initialize property NMRShielding */
void Molecule::initNMRShielding(int k) {
    if (this->nmr.size() == 0) initNuclearProperties(getNNuclei());
    if (this->nmr[k] != nullptr) MSG_ERROR("NMR shielding tensor already initialized");
    this->nmr[k] = std::make_unique<NMRShielding>(getNuclei()[k]);
}

/** @brief Initialize property HyperFineCoupling */
void Molecule::initHyperFineCoupling(int k) {
    if (this->hfcc.size() == 0) initNuclearProperties(getNNuclei());
    if (this->hfcc[k] != nullptr) MSG_ERROR("HyperFine coupling tensor already initialized");
    this->hfcc[k] = std::make_unique<HyperFineCoupling>(getNuclei()[k]);
}

/** @brief Initialize property SpinSpinCoupling */
void Molecule::initSpinSpinCoupling(int k, int l) {
    if (this->sscc.size() == 0) initNuclearProperties(getNNuclei());
    if (this->sscc[k][l] != nullptr) MSG_ERROR("Spin-spin coupling tensor already initialized");
    this->sscc[k][l] = std::make_unique<SpinSpinCoupling>(getNuclei()[k], getNuclei()[l]);
}

/** @brief Initialize property Polarizability */
void Molecule::initPolarizability(double omega) {
    for (auto &i : this->polarizability) {
        const Polarizability &pol_i = *i;
        auto omega_i = pol_i.getFrequency();
        if (std::abs(omega_i - omega) < mrcpp::MachineZero) {
            MSG_ERROR("Polarizability already initialized");
            return;
        }
    }
    this->polarizability.push_back(std::make_unique<Polarizability>(omega));
}

/** @brief Initialize property OpticalRotation */
void Molecule::initOpticalRotation(double omega) {
    NOT_IMPLEMENTED_ABORT;
    /*
    for (auto i = 0; i < this->optical_rotation.size(); i++) {
        OpticalRotation &optrot_i = *this->optical_rotation[i];
        auto omega_i = optrot_i.getFrequency();
        if (fabs(omega_i - omega) < MachineZero) {
            MSG_ERROR("OpticalRotation already initialized");
            return;
        }
    }
    OpticalRotation *optrot = new OpticalRotation(omega);
    this->optical_rotation.push_back(optrot);
    */
}

/** @brief Return property SCFEnergy */
SCFEnergy &Molecule::getSCFEnergy() {
    if (this->energy == nullptr) this->energy = std::make_unique<SCFEnergy>();
    return *this->energy;
}

/** @brief Return property DipoleMoment */
DipoleMoment &Molecule::getDipoleMoment() {
    if (this->dipole == nullptr) MSG_ERROR("Uninitialized dipole moment");
    return *this->dipole;
}

/** @brief Return property QuadrupoleMoment */
QuadrupoleMoment &Molecule::getQuadrupoleMoment() {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (this->quadrupole == nullptr) MSG_ERROR("Uninitialized quadrupole moment");
    return *this->quadrupole;
    */
}

/** @brief Return property GeometryDerivatives */
GeometryDerivatives &Molecule::getGeometryDerivatives() {
    if (this->geomderiv == nullptr) MSG_ERROR("Uninitialized geometry derivatives");
    return *this->geomderiv;
}

/** @brief Return property Magnetizability */
Magnetizability &Molecule::getMagnetizability() {
    if (this->magnetizability == nullptr) MSG_ERROR("Uninitialized magnetizability");
    return *this->magnetizability;
}

/** @brief Return property NMRShielding */
NMRShielding &Molecule::getNMRShielding(int k) {
    if (this->nmr[k] == nullptr) MSG_ERROR("Uninitialized NMR shielding tensor " << k);
    return *this->nmr[k];
}

/** @brief Return property HyperFineCoupling */
HyperFineCoupling &Molecule::getHyperFineCoupling(int k) {
    if (this->hfcc[k] == nullptr) MSG_ERROR("Uninitialized hyperfine coupling tensor " << k);
    return *this->hfcc[k];
}

/** @brief Return property SpinSpinCoupling */
SpinSpinCoupling &Molecule::getSpinSpinCoupling(int k, int l) {
    if (this->sscc[k][l] == nullptr) MSG_ERROR("Uninitialized spin-spin coupling tensor " << k << " " << l);
    return *this->sscc[k][l];
}

/** @brief Return property Polarizability */
Polarizability &Molecule::getPolarizability(double omega) {
    auto idx = -1;
    for (auto i = 0; i < this->polarizability.size(); i++) {
        auto omega_i = this->polarizability[i]->getFrequency();
        if (std::abs(omega_i - omega) < mrcpp::MachineZero) {
            idx = i;
            break;
        }
    }
    if (idx < 0) MSG_ERROR("Uninitialized polarizability");
    return *(this->polarizability[idx]);
}

/** @brief Return property OpticalRotation */
OpticalRotation &Molecule::getOpticalRotation(double omega) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Return number of electrons */
int Molecule::getNElectrons() const {
    auto totZ = 0;
    for (auto i = 0; i < getNNuclei(); i++) totZ += getNuclei()[i].getElement().getZ();
    return totZ - this->charge;
}

/** @brief Compute nuclear center of mass */
Coord<3> Molecule::calcCenterOfMass() const {
    Coord<3> COM;
    COM.fill(0.0);

    auto M = 0.0;
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        const auto &r_i = nuc.getCoord();
        const auto m_i = nuc.getElement().getMass();
        for (auto d = 0; d < 3; d++) COM[d] += r_i[d] * m_i;
        M += m_i;
    }
    for (auto d = 0; d < 3; d++) COM[d] *= 1.0 / M;
    return COM;
}

/** @brief Compute nuclear center of charge */
Coord<3> Molecule::calcCenterOfCharge() const {
    Coord<3> COC;
    COC.fill(0.0);

    auto Z = 0.0;
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        const auto &r_i = nuc.getCoord();
        const auto z_i = nuc.getElement().getZ();
        for (auto d = 0; d < 3; d++) COC[d] += r_i[d] * z_i;
        Z += z_i;
    }
    for (auto d = 0; d < 3; d++) COC[d] *= 1.0 / Z;
    return COC;
}

/** @brief Read nuclear coordinates from xyz file
 *
 * First entry in file is number of atoms:
 *
 * nAtoms
 * symbol   x_coord     y_coord     z_coord
 * symbol   x_coord     y_coord     z_coord
 * symbol   x_coord     y_coord     z_coord
 *
 */
void Molecule::readCoordinateFile(const std::string &coord_file) {
    std::ifstream ifs(coord_file.c_str());
    if (not ifs) MSG_FATAL("Failed to open coordinate file: " << coord_file);

    int nNuclei;
    Coord<3> coord;
    std::string sym;
    ifs >> nNuclei;
    for (auto i = 0; i < nNuclei; i++) {
        ifs >> sym;
        ifs >> coord[0];
        ifs >> coord[1];
        ifs >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
    ifs.close();
}

/** @brief Read nuclear coordinates from vector of strings
 *
 * Each entry in the vector of strings contains one atom:
 *
 *      "symbol   x_coord     y_coord     z_coord"
 *
 */
void Molecule::readCoordinateString(const std::vector<std::string> &coord_str) {
    Coord<3> coord;
    std::string sym;
    for (const auto &i : coord_str) {
        std::stringstream ss;
        ss.str(i);
        ss >> sym;
        ss >> coord[0];
        ss >> coord[1];
        ss >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
}

/** @brief Pretty output of molecular geometry */
void Molecule::printGeometry() const {
    Printer::printHeader(0, "Molecule");
    println(0, " Nr  Element             x             y             z      ");
    Printer::printSeparator(0, '-');
    auto oldPrec = Printer::setPrecision(5);

    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        const auto &coord = nuc.getCoord();
        std::stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, std::setw(3) << i + 1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, std::setw(21) << coord[0]);
        printout(0, std::setw(14) << coord[1]);
        printout(0, std::setw(14) << coord[2] << std::endl);
    }
    Printer::printSeparator(0, '-');
    printout(0, " Center of mass: ");
    Coord<3> COM = calcCenterOfMass();
    printout(0, std::setw(14) << COM[0]);
    printout(0, std::setw(14) << COM[1]);
    printout(0, std::setw(14) << COM[2] << std::endl);
    Printer::setPrecision(oldPrec);
    Printer::printSeparator(0, '=', 2);
}

/** @brief Pretty output of molecular properties
 *
 * Only properties that have been initialized will be printed.
 */
void Molecule::printProperties() const {
    const auto &Phi = *getOrbitals();
    const auto &F_mat = getFockMatrix();
    orbital::print_eigenvalues(Phi, F_mat);

    if (this->energy != nullptr) println(0, *this->energy);
    if (this->dipole != nullptr) println(0, *this->dipole);
    if (this->geomderiv != nullptr) println(0, *this->geomderiv);
    if (this->magnetizability != nullptr) println(0, *this->magnetizability);
    for (auto &pol : this->polarizability) {
        if (pol != nullptr) println(0, *pol);
    }
    for (auto &nmr_k : this->nmr) {
        if (nmr_k != nullptr) println(0, *nmr_k);
    }
    for (auto &hfcc_k : this->hfcc) {
        if (hfcc_k != nullptr) println(0, *hfcc_k);
    }
    for (auto &sscc_k : this->sscc) {
        for (auto &sscc_kl : sscc_k) {
            if (sscc_kl != nullptr) println(0, *sscc_kl);
        }
    }
}

} // namespace mrchem
